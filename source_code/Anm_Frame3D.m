%% Anm_Frame3D Class
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <anm.html *Anm*> in the <main.html LESM> (Linear Elements
% Structure Model) program. This sub-class implements abstract methods,
% declared in super-class *Anm*, that deal with 3D frame analysis model
% of linear structure elements.
%
%% Frame 3D Analysis
% These are the basic assumptions of a 2D frame model:
%%%
% * Frame elements are usually rigidly connected at joints.
%   However, a frame element might have a hinge (rotation liberation)
%   at an end or hinges at both ends.
% * It is assumed that a hinge in a 3D frame element releases continuity
%   of rotation in all directions.
% * Internal forces at any cross-section of a 3D frame element are:
%   axial force, shear forces (local y direction and local z direction),
%   bending moments (about local y direction and local z direction),
%   and torsion moment (about x direction).
% * Each node of a 3D frame model has six d.o.f.'s: displacements in
%   X, Y and Z directions, and rotations about X, Y and Z directions.
%
%% Class definition
% Definition of sub-class *Anm_Frame3D* derived from super-class <anm.html Anm>.
classdef Anm_Frame3D < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Frame3D()
            include_constants;
            anm = anm@Anm(FRAME3D_ANALYSIS,6);
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
            % rot = [ T 0 0 0
            %         0 T 0 0
            %         0 0 T 0
            %         0 0 0 T ]
            rot = blkdiag(T,T,T,T);
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
                
                % Check for fixed translation in global Z direction
                if drv.nodes(n).ebc(3) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(3,n) = 1;
                end
                
                % Check for fixed rotation about global X direction
                if drv.nodes(n).ebc(4) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(4,n) = 1;
                end
                
                % Check for fixed rotation about global Y direction
                if drv.nodes(n).ebc(5) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(5,n) = 1;
                end
                
                % Check for fixed rotation about global Z direction
                if drv.nodes(n).ebc(6) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(6,n) = 1;
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
                    
                    % Add prescribed displacement in global Z direction
                    id = drv.ID(3,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(3) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(3);
                    end
                    
                    % Add prescribed rotation about global X direction
                    id = drv.ID(4,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(4) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(4);
                    end
                    
                    % Add prescribed rotation about global Y direction
                    id = drv.ID(5,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(5) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(5);
                    end
                    
                    % Add prescribed rotation about global Z direction
                    id = drv.ID(6,n);
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
            
            % Compute torsion stiffness coefficients
            ket = elem.torsionStiffCoeff();
            
            % Compute flexural stiffness coefficients
            kef_XY = elem.flexuralStiffCoeff_XY();
            kef_XZ = elem.flexuralStiffCoeff_XZ();
            
            % Assemble element stiffness matrix in local system
            kel = [ kea(1,1)  0            0           0         0            0           kea(1,2)  0            0           0         0            0;
                    0         kef_XY(1,1)  0           0         0            kef_XY(1,2) 0         kef_XY(1,3)  0           0         0            kef_XY(1,4);
                    0         0            kef_XZ(1,1) 0         kef_XZ(1,2)  0           0         0            kef_XZ(1,3) 0         kef_XZ(1,4)  0;
                    0         0            0           ket(1,1)  0            0           0         0            0           ket(1,2)  0            0;
                    0         0            kef_XZ(2,1) 0         kef_XZ(2,2)  0           0         0            kef_XZ(2,3) 0         kef_XZ(2,4)  0;
                    0         kef_XY(2,1)  0           0         0            kef_XY(2,2) 0         kef_XY(2,3)  0           0         0            kef_XY(2,4);
                    kea(2,1)  0            0           0         0            0           kea(2,2)  0            0           0         0            0;
                    0         kef_XY(3,1)  0           0         0            kef_XY(3,2) 0         kef_XY(3,3)  0           0         0            kef_XY(3,4);
                    0         0            kef_XZ(3,1) 0         kef_XZ(3,2)  0           0         0            kef_XZ(3,3) 0         kef_XZ(3,4)  0;
                    0         0            0           ket(2,1)  0            0           0         0            0           ket(2,2)  0            0;
                    0         0            kef_XZ(4,1) 0         kef_XZ(4,2)  0           0         0            kef_XZ(4,3) 0         kef_XZ(4,4)  0;
                    0         kef_XY(4,1)  0           0         0            kef_XY(4,2) 0         kef_XY(4,3)  0           0         0            kef_XY(4,4) ];
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
                    
                    % Add applied force in global Z direction
                    id = drv.ID(3,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(3);
                    
                    % Add applied moment about global X direction
                    id = drv.ID(4,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(4);
                    
                    % Add applied moment about global Y direction
                    id = drv.ID(5,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(5);
                    
                    % Add applied moment about global Z direction
                    id = drv.ID(6,n);
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
            fef_XY = load.flexuralDistribLoadFEF_XY();
            fef_XZ = load.flexuralDistribLoadFEF_XZ();
            
            % Assemble element fixed end force (FEF) vector in local system
            fel = [ fea(1);
                    fef_XY(1);
                    fef_XZ(1);
                    0;
                    fef_XZ(2);
                    fef_XY(2);
                    fea(2);
                    fef_XY(3);
                    fef_XZ(3);
                    0;
                    fef_XZ(4);
                    fef_XY(4)];
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
            fef_XY = load.flexuralThermalLoadFEF_XY();
            fef_XZ = load.flexuralThermalLoadFEF_XZ();
            
            % Assemble element fixed end force (FEF) vector in local system
            fel = [ fea(1);
                    fef_XY(1);
                    fef_XZ(1);
                    0;
                    fef_XZ(2);
                    fef_XY(2);
                    fea(2);
                    fef_XY(3);
                    fef_XZ(3);
                    0;
                    fef_XZ(4);
                    fef_XY(4)];
        end
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        %  axial_force(1,2)
        %   Ni = axial_force(1) - init value
        %   Nf = axial_force(2) - final value
        %  shear_force(1,2)
        %   Qi = shear_force(1) - init value
        %   Qf = shear_force(2) - final value
        %  torsion_moment(1,2)
        %   Ti = torsion_moment(1) - init value
        %   Tf = torsion_moment(2) - final value
        %  bending_moment(1,2)
        %   Mi = bending_moment(1) - init value
        %   Mf = bending_moment(2) - final value
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function initIntForce(~,elem)
            elem.axial_force      = zeros(1,2);
            elem.shear_force_Y    = zeros(1,2);
            elem.shear_force_Z    = zeros(1,2);
            elem.torsion_moment   = zeros(1,2);
            elem.bending_moment_Y = zeros(1,2);
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
            elem.axial_force(2)      = elem.axial_force(2)      + fel(7);
            elem.shear_force_Y(1)    = elem.shear_force_Y(1)    + fel(2);
            elem.shear_force_Y(2)    = elem.shear_force_Y(2)    + fel(8);
            elem.shear_force_Z(1)    = elem.shear_force_Z(1)    + fel(3);
            elem.shear_force_Z(2)    = elem.shear_force_Z(2)    + fel(9);
            elem.torsion_moment(1)   = elem.torsion_moment(1)   + fel(4);
            elem.torsion_moment(2)   = elem.torsion_moment(2)   + fel(10);
            elem.bending_moment_Y(1) = elem.bending_moment_Y(1) + fel(5);
            elem.bending_moment_Y(2) = elem.bending_moment_Y(2) + fel(11);
            elem.bending_moment_Z(1) = elem.bending_moment_Z(1) + fel(6);
            elem.bending_moment_Z(2) = elem.bending_moment_Z(2) + fel(12);
        end
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        %  intDispl(1,:) -> du (axial displacement)
        %  intDispl(2,:) -> dv (transversal displacement in local y-axis)
        %  intDispl(3,:) -> dw (transversal displacement in local z-axis)
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function initIntDispl(~,elem)
            elem.intDispl = zeros(3,50);
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
            Nw = elem.flexuralDisplShapeFcnVector_XZ(x);
            
            % Assemble displacement shape function matrix
            N = [ Nu(1)  0      0      0      0      0      Nu(2)  0      0      0      0      0;
                  0      Nv(1)  0      0      0      Nv(2)  0      Nv(3)  0      0      0      Nv(4);
                  0      0      Nw(1)  0      Nw(2)  0      0      0      Nw(3)  0      Nw(4)  0 ];
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements vector in local system,
        % in a given cross-section position, for the local analysis from
        % element loads (distributed loads and thermal loads).
        % Output:
        %  del: a 3x1 vector of with element internal displacements:
        %      du -> axial displacement
        %      dv -> transversal displacement in local y-axis direction
        %      dw -> transversal displacement in local z-axis direction
        %      x  -> cross-section position on element local x-axis
        %  del = [ du(x);
        %          dv(x);
        %          dw(x) ]
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function del = lclAnlIntDispl(~,elem,x)
            % Initialize displacements vector resulting from local analysis
            del = [ 0;
                    0;
                    0 ];
            
            % Add the contribution of axial and transversal displacements
            % resulting from distributed loads
            if (~isempty(elem.load.uniformLcl)) || (~isempty(elem.load.linearLcl))
                del(1) = elem.load.axialDistribLoadDispl(x);
                del(2) = elem.load.flexuralDistribLoadDispl_XY(x);
                del(3) = elem.load.flexuralDistribLoadDispl_XZ(x);
            end
            
            % Add the contribution of axial and transversal displacements
            % resulting from thermal loads
            if (elem.load.tempVar_X ~= 0) || (elem.load.tempVar_Y ~= 0) || (elem.load.tempVar_Z ~= 0)
                del(1) = del(1) + elem.load.axialThermalLoadDispl(x);
                del(2) = del(2) + elem.load.flexuralThermalLoadDispl_XY(x);
                del(3) = del(3) + elem.load.flexuralThermalLoadDispl_XZ(x);
            end
        end
    end
end