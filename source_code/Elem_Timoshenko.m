%% Elem_Timoshenko Class
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <elem.html *Elem*> in the <main.html LESM> (Linear Elements
% Structure Model) program. This sub-class implements abstract methods,
% declared in super-class *Elem*, that deal with Timoshenko flexural
% behavior of linear elements.
%
%% Timoshenko Beam Theory
%
% In Timoshenko flexural behavior, shear deformation is considered in an
% approximated manner. Bending of a linear structure element is such that its
% cross-section remains plane but it is not normal to the element longitudinal
% axis.
%
%% Class definition
% Definition of sub-class *Elem_Timoshenko* derived from super-class <elem.html Elem>.
classdef Elem_Timoshenko < Elem
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem_Timoshenko(anm,mat,sec,nodes,hi,hf,vz,gle,af,sfy,sfz,bmy,bmz,tm,intdispl)
            include_constants;
            elem = elem@Elem(MEMBER_TIMOSHENKO);
            
            if (nargin > 0)
                % Get nodal coordinates
                xi = nodes(1).coord(1);
                yi = nodes(1).coord(2);
                zi = nodes(1).coord(3);
                xf = nodes(2).coord(1);
                yf = nodes(2).coord(2);
                zf = nodes(2).coord(3);
                
                % Calculate element length
                dx = xf - xi;
                dy = yf - yi;
                dz = zf - zi;
                L  = sqrt(dx^2 + dy^2 + dz^2);
                
                % Calculate element cosines with global axes
                cx = dx/L;
                cy = dy/L;
                cz = dz/L;
                
                % Create load object
                load = Lelem_Timoshenko(elem);

                % Set properties
                elem.anm = anm;
                elem.load = load;
                elem.material = mat;
                elem.section = sec;
                elem.nodes = nodes;
                elem.hingei = hi;
                elem.hingef = hf;
                elem.length = L;
                elem.cosine_X = cx;
                elem.cosine_Y = cy;
                elem.cosine_Z = cz;
                elem.vz = vz;
                elem.gle = gle;
                elem.axial_force = af;
                elem.shear_force_Y = sfy;
                elem.shear_force_Z = sfz;
                elem.bending_moment_Y = bmy;
                elem.bending_moment_Z = bmz;
                elem.torsion_moment = tm;
                elem.intDispl = intdispl;
                
                % Compute and set basis rotation transformation matrix and
                % d.o.f. rotation transformation matrix
                elem.T = elem.rotTransMtx();
                elem.rot = anm.gblToLocElemRotMtx(elem);
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <elem.html *Elem*>.
    methods
        %------------------------------------------------------------------
        % Generates element flexural stiffness coefficient matrix in local
        % xy-plane.
        % Output:
        %  kef: a 4x4 matrix with flexural stiffness coefficients:
        %      D -> transversal displacement in local y-axis
        %      R -> rotation about local z-axis
        %      ini -> initial node
        %      end -> end node
        %   kef = [ k_Dini_Dini  k_Dini_Rini  k_Dini_Dend  k_Dini_Rend;
        %           k_Rini_Dini  k_Rini_Rini  k_Rini_Dend  k_Dini_Rend;
        %           k_Dend_Dini  k_Dend_Rini  k_Dend_Dend  k_Dend_Rend;
        %           k_Rend_Dini  k_Rend_Rini  k_Rend_Dend  k_Rend_Rend ]
        function kef = flexuralStiffCoeff_XY(elem)
            include_constants;
            
            E  = elem.material.elasticity;
            G  = elem.material.shear;
            As = elem.section.area_y;
            I  = elem.section.inertia_z;
            L  = elem.length;
            
            L2  = L^2;
            GAs = G * As;
            EI  = E * I;
            
            % Timoshenko parameters
            Omega  = EI / (GAs * L2);
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            laL  = lambda * L;
            laL2 = laL    * L;
            laL3 = laL2   * L;
            
            muL  = mu   * L;
            muL2 = muL  * L;
            muL3 = muL2 * L;
            
            if (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [ 12*EI/muL3   6*EI/muL2        -12*EI/muL3   6*EI/muL2;
                         6*EI/muL2   4*lambda*EI/muL   -6*EI/muL2   2*gamma*EI/muL;
                       -12*EI/muL3  -6*EI/muL2         12*EI/muL3  -6*EI/muL2;
                         6*EI/muL2   2*gamma*EI/muL    -6*EI/muL2   4*lambda*EI/muL ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [  3*EI/laL3   0                 -3*EI/laL3   3*EI/laL2;
                         0           0                  0           0;
                        -3*EI/laL3   0                  3*EI/laL3  -3*EI/laL2;
                         3*EI/laL2   0                 -3*EI/laL2   3*EI/laL ];
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                
                kef = [  3*EI/laL3   3*EI/laL2         -3*EI/laL3   0;
                         3*EI/laL2   3*EI/laL          -3*EI/laL2   0;
                        -3*EI/laL3  -3*EI/laL2          3*EI/laL3   0;
                         0           0                  0           0 ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                
                kef = [  0           0                  0           0;
                         0           0                  0           0;
                         0           0                  0           0;
                         0           0                  0           0 ];
                
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural stiffness coefficient matrix in local
        % xz-plane.
        % Output:
        %  kef: a 4x4 matrix with flexural stiffness coefficients:
        %      D -> transversal displacement in local z-axis
        %      R -> rotation about local y-axis
        %      ini -> initial node
        %      end -> end node
        %   kef = [ k_Dini_Dini  k_Dini_Rini  k_Dini_Dend  k_Dini_Rend;
        %           k_Rini_Dini  k_Rini_Rini  k_Rini_Dend  k_Dini_Rend;
        %           k_Dend_Dini  k_Dend_Rini  k_Dend_Dend  k_Dend_Rend;
        %           k_Rend_Dini  k_Rend_Rini  k_Rend_Dend  k_Rend_Rend ]
        function kef = flexuralStiffCoeff_XZ(elem)
            include_constants;
            
            E  = elem.material.elasticity;
            G  = elem.material.shear;
            As = elem.section.area_z;
            I  = elem.section.inertia_y;
            L  = elem.length;
            
            L2  = L^2;
            GAs = G * As;
            EI  = E * I;
            
            % Timoshenko parameters
            Omega  = EI / (GAs * L2); 
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            laL  = lambda * L;
            laL2 = laL    * L;
            laL3 = laL2   * L;
            
            muL  = mu   * L;
            muL2 = muL  * L;
            muL3 = muL2 * L;
            
            if (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [ 12*EI/muL3  -6*EI/muL2        -12*EI/muL3  -6*EI/muL2;
                        -6*EI/muL2   4*lambda*EI/muL    6*EI/muL2   2*gamma*EI/muL;
                       -12*EI/muL3   6*EI/muL2         12*EI/muL3   6*EI/muL2;
                        -6*EI/muL2   2*gamma*EI/muL     6*EI/muL2   4*lambda*EI/muL ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [  3*EI/laL3   0                 -3*EI/laL3  -3*EI/laL2;
                         0           0                  0           0;
                        -3*EI/laL3   0                  3*EI/laL3   3*EI/laL2;
                        -3*EI/laL2   0                  3*EI/laL2   3*EI/laL ];
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                
                kef = [  3*EI/laL3  -3*EI/laL2         -3*EI/laL3   0;
                        -3*EI/laL2   3*EI/laL           3*EI/laL2   0;
                        -3*EI/laL3   3*EI/laL2          3*EI/laL3   0;
                         0           0                  0           0 ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                
                kef = [  0           0                  0           0;
                         0           0                  0           0;
                         0           0                  0           0;
                         0           0                  0           0 ];
                
            end
        end
        
        %------------------------------------------------------------------
        % Evaluates element flexural displacement shape functions in local
        % xy-plane at a given cross-section position.
        % Flexural shape functions give the transversal displacement of an
        % element subjected to unit nodal displacements and rotations.
        % Output:
        %  Nv: a 1x4 vector with flexural displacement shape functions:
        %      D -> transversal displacement in local y-axis
        %      R -> rotation about local z-axis
        %      ini -> initial node
        %      end -> end node
        %   Nv = [ N_Dini  N_Rini  N_Dend  N_Rend ]
        % Input arguments:
        %  x: element cross-section position in local x-axis
        function Nv = flexuralDisplShapeFcnVector_XY(elem,x)
            include_constants;
            
            E  = elem.material.elasticity;
            G  = elem.material.shear;
            As = elem.section.area_y;
            I  = elem.section.inertia_z;
            L  = elem.length;
            
            L2  = L^2;
            L3  = L2 * L;
            x2  = x^2;
            x3  = x * x2;
            GAs = G * As;
            EI  = E * I;
            
            % Timoshenko parameters
            Omega  = EI / (GAs * L2);
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            if  (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                Nv1 = 1 - 12*Omega*x/(L*mu) - 3*x2/(L2*mu) + 2*x3/(L3*mu);
                Nv2 = x - 6*Omega*x/mu - 2*lambda*x2/(L*mu) + x3/(L2*mu);
                Nv3 = 12*Omega*x/(L*mu) + 3*x2/(L2*mu) - 2*x3/(L3*mu);
                Nv4 = -6*Omega*x/mu - gamma*x2/(L*mu) + x3/(L2*mu);
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                Nv1 = 1 - 3*x/(2*L*lambda) - 3*Omega*x/(L*lambda) + x3/(2*L3*lambda);
                Nv2 = 0;
                Nv3 = 3*x/(2*L*lambda) + 3*Omega*x/(L*lambda) - x3/(2*L3*lambda);
                Nv4 = -gamma*x/(2*lambda) - 3*Omega*x/lambda + x3/(2*L2*lambda);
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                Nv1 = 1 - 3*Omega*x/(L*lambda) - 3*x2/(2*L2*lambda) + x3/(2*L3*lambda);
                Nv2 = x - 3*Omega*x/lambda - 3*x2/(2*L*lambda) + x3/(2*L2*lambda);
                Nv3 = 3*Omega*x/(L*lambda) + 3*x2/(2*L2*lambda) - x3/(2*L3*lambda);
                Nv4 = 0;
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                Nv1 = 1 - x/L;
                Nv2 = 0;
                Nv3 = x/L;
                Nv4 = 0;
            end
            
            Nv = [ Nv1  Nv2  Nv3  Nv4 ];
        end
        
        %------------------------------------------------------------------
        % Evaluates element flexural displacement shape functions in local
        % xz-plane at a given cross-section position.
        % Flexural shape functions give the transversal displacement of an
        % element subjected to unit nodal displacements and rotations.
        % Output:
        %  Nw: a 1x4 vector with flexural displacement shape functions:
        %      D -> transversal displacement in local z-axis
        %      R -> rotation about local y-axis
        %      ini -> initial node
        %      end -> end node
        %   Nw = [ N_Dini  N_Rini  N_Dend  N_Rend ]
        % Input arguments:
        %  x: element cross-section position in local x-axis
        function Nw = flexuralDisplShapeFcnVector_XZ(elem,x)
            include_constants;
            
            E  = elem.material.elasticity;
            G  = elem.material.shear;
            As = elem.section.area_z;
            I  = elem.section.inertia_y;
            L  = elem.length;
            
            L2  = L^2;
            L3  = L2 * L;
            x2  = x^2;
            x3  = x * x2;
            GAs = G * As;
            EI  = E * I;
            
            % Timoshenko parameters
            Omega  = EI / (GAs * L2);
            lambda = 1 + 3  * Omega;
            mu     = 1 + 12 * Omega;
            gamma  = 1 - 6  * Omega;
            
            if  (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                Nw1 = 1 - 12*Omega*x/(L*mu) - 3*x2/(L2*mu) + 2*x3/(L3*mu);
                Nw2 = -x + 6*Omega*x/mu + 2*lambda*x2/(L*mu) - x3/(L2*mu);
                Nw3 = 12*Omega*x/(L*mu) + 3*x2/(L2*mu) - 2*x3/(L3*mu);
                Nw4 = 6*Omega*x/mu + gamma*x2/(L*mu) - x3/(L2*mu);
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                Nw1 = 1 - 3*x/(2*L*lambda) - 3*Omega*x/(L*lambda) + x3/(2*L3*lambda);
                Nw2 = 0;
                Nw3 = 3*x/(2*L*lambda) + 3*Omega*x/(L*lambda) - x3/(2*L3*lambda);
                Nw4 = gamma*x/(2*lambda) + 3*Omega*x/lambda - x3/(2*L2*lambda);
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                Nw1 = 1 - 3*Omega*x/(L*lambda) - 3*x2/(2*L2*lambda) + x3/(2*L3*lambda);
                Nw2 = -x + 3*Omega*x/lambda + 3*x2/(2*L*lambda) - x3/(2*L2*lambda);
                Nw3 = 3*Omega*x/(L*lambda) + 3*x2/(2*L2*lambda) - x3/(2*L3*lambda);
                Nw4 = 0;
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                Nw1 = 1 - x/L;
                Nw2 = 0;
                Nw3 = x/L;
                Nw4 = 0;
            end
            
            Nw = [ Nw1  Nw2  Nw3  Nw4 ];
        end
    end
end