%% Anm_Frame2D Class
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <anm.html *Anm*> in the <main.html LESM> (Linear Elements
% Structure Model) program. This sub-class implements abstract methods,
% declared in super-class *Anm*, that deal with 2D frame analysis model
% of linear structure elements.
%
%% Frame 2D Analysis
% These are the basic assumptions of a 2D frame model:
%%%
% * Frame elements are usually rigidly connected at joints.
%   However, a frame element might have a hinge (rotation liberation)
%   at an end or hinges at both ends.
% * It is assumed that a hinge in a 2D frame element releases continuity of
%   rotation in the out-of-plane direction.
% * A 2D frame model is considered to be laid in the XY-plane, with only
%   in-plane behavior, that is, there is no displacement transversal to
%   the frame plane.
% * Internal forces at any cross-section of a 2D frame element are:
%   axial force (local x direction), shear force (local y direction),
%   and bending moment (about local z direction).
% * Each node of a 2D frame model has three d.o.f.'s: a horizontal
%   displacement in X direction, a vertical displacement in Y
%   direction, and a rotation about the Z-axis.
%
%% Class definition
% Definition of sub-class *Anm_Frame2D* derived from super-class <anm.html Anm>.
classdef Anm_Frame2D < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Frame2D()
            include_constants;
            anm = anm@Anm(FRAME2D_ANALYSIS,3);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <anm.html *Anm*>.
    methods
        %------------------------------------------------------------------
        % Assembles element d.o.f. (degree of freedom) rotation transformation
        % matrix from global system to local system.
        % Output:
        %  rot: rotation transformation matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function rot = gblToLocElemRotMtx(~,elem)
            % Get 3x3 basis rotation transformation matrix
            T = elem.T;
            
            % Assemble element d.o.f. rotation transformation matrix
            % rot = [ T 0
            %         0 T ]
            rot = blkdiag(T,T);
        end
        
        %------------------------------------------------------------------
        % Initializes global d.o.f (degree of freedom) numbering ID matrix
        % with ones and zeros, and counts total number of equations of free
        % d.o.f.'s and total number of equations of fixed d.o.f.'s.
        %  ID matrix initialization:
        %  if ID(k,n) = 0, d.o.f. k of node n is free.
        %  if ID(k,n) = 1, d.o.f. k of node n is constrained by support.
        % Input arguments:
        %  drv: handle to an object of the Drv class
        function setupDOFNum(anm,drv)
            include_constants;
            
            % Dimension global d.o.f. numbering ID matrix
            drv.ID = zeros(anm.ndof,drv.nnp);
            
            % Initialize number of fixed d.o.f.
            drv.neqfixed = 0;
            
            % Count number of fixed d.o.f. and setup ID matrix
            for n = 1:drv.nnp
                
                % Check for fixed translation in global X direction
                if drv.nodes(n).ebc(1) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(1,n) = 1;
                end
                
                % Check for fixed translation in global Y direction
                if drv.nodes(n).ebc(2) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(2,n) = 1;
                end
                
                % Check for fixed rotation about global Z direction
                if drv.nodes(n).ebc(6) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(3,n) = 1;
                end
            end
            
            % Compute total number of free d.o.f.
            drv.neqfree = drv.neq - drv.neqfixed;
        end
        
        %------------------------------------------------------------------
        % Adds prescribed displacements (known support settlement values)
        % to global displacement vector.
        % Avoids storing a prescribed displacement component in a position
        % of the global displacement vector that corresponds to a free d.o.f.
        % Input arguments:
        %  drv: handle to an object of the Drv class
        function setupPrescDispl(~,drv)
            for n = 1:drv.nnp
                if ~isempty(drv.nodes(n).prescDispl)
                    
                    % Add prescribed displacement in global X direction
                    id = drv.ID(1,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(1) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(1);
                    end
                    
                    % Add prescribed displacement in global Y direction
                    id = drv.ID(2,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(2) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(2);
                    end
                    
                    % Add prescribed rotation about global Z direction
                    id = drv.ID(3,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(6) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(6);
                    end
                    
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element stiffness matrix in local system.
        % Output:
        %  kel: target element stiffness matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function kel = elemLocStiffMtx(~,elem)
            % Compute axial stiffness coefficients
            kea = elem.axialStiffCoeff();
            
            % Compute flexural stiffness coefficients
            kef = elem.flexuralStiffCoeff_XY();
            
            % Assemble element stiffness matrix in local system
            kel = [ kea(1,1)  0         0         kea(1,2)  0         0;
                    0         kef(1,1)  kef(1,2)  0         kef(1,3)  kef(1,4);
                    0         kef(2,1)  kef(2,2)  0         kef(2,3)  kef(2,4);
                    kea(2,1)  0         0         kea(2,2)  0         0;
                    0         kef(3,1)  kef(3,2)  0         kef(3,3)  kef(3,4);
                    0         kef(4,1)  kef(4,2)  0         kef(4,3)  kef(4,4) ];
        end
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  drv: handle to an object of the Drv class
        function nodalLoads(~,drv)
            for n = 1:drv.nnp
                if ~isempty(drv.nodes(n).nodalLoad)
                    
                    % Add applied force in global X direction
                    id = drv.ID(1,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(1);
                    
                    % Add applied force in global Y direction
                    id = drv.ID(2,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(2);
                    
                    % Add applied moment about global Z direction
                    id = drv.ID(3,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(6);
                    
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied distributed load.
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        function fel = elemLocDistribLoadFEF(~,load)
            % Compute axial fixed end force components
            fea = load.axialDistribLoadFEF();
            
            % Compute flexural (transversal) fixed end force components
            fef = load.flexuralDistribLoadFEF_XY();
            
            % Assemble element fixed end force (FEF) vector in local system
            fel = [ fea(1);
                    fef(1);
                    fef(2);
                    fea(2);
                    fef(3);
                    fef(4) ];
        end
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied thermal load (temperature variation).
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        function fel = elemLocThermalLoadFEF(~,load)
            % Compute axial fixed end force components
            fea = load.axialThermalLoadFEF();
            
            % Compute flexural (transversal) fixed end force components
            fef = load.flexuralThermalLoadFEF_XY();
            
            % Assemble element fixed end force (FEF) vector in local system
            fel = [ fea(1);
                    fef(1);
                    fef(2);
                    fea(2);
                    fef(3);
                    fef(4) ];
        end
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        %  axial_force(1,2)
        %   Ni = axial_force(1) - init value
        %   Nf = axial_force(2) - final value
        %  shear_force(1,2)
        %   Qi = shear_force(1) - init value
        %   Qf = shear_force(2) - final value
        %  bending_moment(1,2)
        %   Mi = bending_moment(1) - init value
        %   Mf = bending_moment(2) - final value
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function initIntForce(~,elem)
            elem.axial_force      = zeros(1,2);
            elem.shear_force_Y    = zeros(1,2);
            elem.bending_moment_Z = zeros(1,2);
        end
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        function assembleIntForce(~,elem,fel)
            elem.axial_force(1)      = elem.axial_force(1)      + fel(1);
            elem.axial_force(2)      = elem.axial_force(2)      + fel(4);
            elem.shear_force_Y(1)    = elem.shear_force_Y(1)    + fel(2);
            elem.shear_force_Y(2)    = elem.shear_force_Y(2)    + fel(5);
            elem.bending_moment_Z(1) = elem.bending_moment_Z(1) + fel(3);
            elem.bending_moment_Z(2) = elem.bending_moment_Z(2) + fel(6);
        end
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        %  intDispl(1,:) -> du (axial displacement)
        %  intDispl(2,:) -> dv (transversal displacement in local y-axis)
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function initIntDispl(~,elem)
            elem.intDispl = zeros(2,50);
        end
        
        %------------------------------------------------------------------
        % Assembles displacement shape function matrix evaluated at a 
        % given cross-section position.
        % Output:
        %  N: displacement shape function matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function N = displShapeFcnMtx(~,elem,x)
            % Compute axial displacement shape functions vector
            Nu = elem.axialDisplShapeFcnVector(x);
            
            % Compute transversal displacement shape functions vector
            Nv = elem.flexuralDisplShapeFcnVector_XY(x);
            
            % Assemble displacement shape function matrix
            N = [ Nu(1)  0      0      Nu(2)  0      0;
                  0      Nv(1)  Nv(2)  0      Nv(3)  Nv(4) ];
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements vector in local system,
        % in a given cross-section position, for the local analysis from
        % element loads (distributed loads and thermal loads).
        % Output:
        %  del: a 2x1 vector of with element internal displacements:
        %      du -> axial displacement
        %      dv -> transversal displacement in local y-axis direction
        %      x  -> cross-section position on element local x-axis
        %  del = [ du(x);
        %          dv(x) ]
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function del = lclAnlIntDispl(~,elem,x)
            % Initialize displacements vector resulting from local analysis
            del = [ 0;
                    0 ];
            
            % Add the contribution of axial and transversal displacements
            % resulting from distributed loads
            if (~isempty(elem.load.uniformLcl)) || (~isempty(elem.load.linearLcl))
                del(1) = elem.load.axialDistribLoadDispl(x);
                del(2) = elem.load.flexuralDistribLoadDispl_XY(x);
            end
            
            % Add the contribution of axial and transversal displacements
            % resulting from thermal loads
            if (elem.load.tempVar_X ~= 0) || (elem.load.tempVar_Y ~= 0)
                del(1) = del(1) + elem.load.axialThermalLoadDispl(x);
                del(2) = del(2) + elem.load.flexuralThermalLoadDispl_XY(x);
            end
        end
    end
end