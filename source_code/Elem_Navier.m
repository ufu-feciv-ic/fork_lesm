%% Elem_Navier Class
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <elem.html *Elem*> in the <main.html LESM> (Linear Elements
% Structure Model) program. This sub-class implements abstract methods,
% declared in super-class *Elem*, that deal with Navier (Euler-Bernoulli)
% flexural behavior of linear elements.
%
%% Navier (Euler-Bernoulli) Beam Theory
%
% In Euler-Bernoulli flexural behavior, it is assumed that there is no
% shear deformation. As consequence, bending of a linear element is such that
% its cross-section remains plane and normal to the element longitudinal axis.
%
%% Class definition
% Definition of sub-class *Elem_Navier* derived from super-class <elem.html Elem>.
classdef Elem_Navier < Elem
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem_Navier(anm,mat,sec,nodes,hi,hf,vz,gle,af,sfy,sfz,bmy,bmz,tm,intdispl)
            include_constants;
            elem = elem@Elem(MEMBER_NAVIER);
            
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
                load = Lelem_Navier(elem);

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
            
            E = elem.material.elasticity;
            I = elem.section.inertia_z;
            L = elem.length;
            
            L2 = L^2;
            L3 = L2 * L;
            EI = E  * I;
            
            if (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [ 12*EI/L3   6*EI/L2  -12*EI/L3   6*EI/L2;
                         6*EI/L2   4*EI/L    -6*EI/L2   2*EI/L;
                       -12*EI/L3  -6*EI/L2   12*EI/L3  -6*EI/L2;
                         6*EI/L2   2*EI/L    -6*EI/L2   4*EI/L ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [  3*EI/L3   0         -3*EI/L3   3*EI/L2;
                         0         0          0         0;
                        -3*EI/L3   0          3*EI/L3  -3*EI/L2;
                         3*EI/L2   0         -3*EI/L2   3*EI/L ];
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                
                kef = [  3*EI/L3   3*EI/L2   -3*EI/L3   0;
                         3*EI/L2   3*EI/L    -3*EI/L2   0;
                        -3*EI/L3  -3*EI/L2    3*EI/L3   0;
                         0         0          0         0 ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)

                kef = [  0         0          0         0;
                         0         0          0         0;
                         0         0          0         0;
                         0         0          0         0 ];
                
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
            
            E = elem.material.elasticity;
            I = elem.section.inertia_y;
            L = elem.length;
            
            L2 = L^2;
            L3 = L2 * L;
            EI = E  * I;
            
            if (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [ 12*EI/L3  -6*EI/L2  -12*EI/L3  -6*EI/L2;
                        -6*EI/L2   4*EI/L     6*EI/L2   2*EI/L;
                       -12*EI/L3   6*EI/L2   12*EI/L3   6*EI/L2;
                        -6*EI/L2   2*EI/L     6*EI/L2   4*EI/L ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                
                kef = [  3*EI/L3   0         -3*EI/L3  -3.*EI/L2;
                         0         0          0         0;
                        -3*EI/L3   0          3*EI/L3   3*EI/L2;
                        -3*EI/L2   0          3*EI/L2   3*EI/L ];
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                
                kef = [  3*EI/L3  -3*EI/L2   -3*EI/L3   0;
                        -3*EI/L2   3*EI/L     3*EI/L2   0;
                        -3*EI/L3   3*EI/L2    3*EI/L3   0;
                         0         0          0         0 ];
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == HINGED_END)
                
                kef = [  0         0          0         0;
                         0         0          0         0;
                         0         0          0         0;
                         0         0          0         0 ];
                
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
            
            L  = elem.length;
            L2 = L  * L;
            L3 = L2 * L;
            x2 = x  * x;
            x3 = x2 * x;
            
            if  (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                Nv1 = 1 - 3*x2/L2 + 2*x3/L3;
                Nv2 = x - 2*x2/L + x3/L2;
                Nv3 = 3*x2/L2 - 2*x3/L3;
                Nv4 = -x2/L + x3/L2;
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                Nv1 = 1 - 3*x/(2*L) + x3/(2*L3);
                Nv2 = 0;
                Nv3 = 3*x/(2*L) - x3/(2*L3);
                Nv4 = -x/2 + x3/(2*L2);
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                Nv1 = 1 - 3*x2/(2*L2) + x3/(2*L3);
                Nv2 = x - 3*x2/(2*L) + x3/(2*L2);
                Nv3 = 3*x2/(2*L2) - x3/(2*L3);
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
            
            L  = elem.length;
            L2 = L  * L;
            L3 = L2 * L;
            x2 = x  * x;
            x3 = x2 * x;
            
            if  (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                Nw1 = 1 - 3*x2/L2 + 2*x3/L3;
                Nw2 = -x + 2*x2/L - x3/L2;
                Nw3 = 3*x2/L2 - 2*x3/L3;
                Nw4 = x2/L - x3/L2;
                
            elseif (elem.hingei == HINGED_END) && (elem.hingef == CONTINUOUS_END)
                Nw1 = 1 - 3*x/(2*L) + x3/(2*L3);
                Nw2 = 0;
                Nw3 = 3*x/(2*L) - x3/(2*L3);
                Nw4 = x/2 - x3/(2*L2);
                
            elseif (elem.hingei == CONTINUOUS_END) && (elem.hingef == HINGED_END)
                Nw1 = 1 - 3*x2/(2*L2) + x3/(2*L3);
                Nw2 = -x + 3*x2/(2*L) - x3/(2*L2);
                Nw3 = 3*x2/(2*L2) - x3/(2*L3);
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