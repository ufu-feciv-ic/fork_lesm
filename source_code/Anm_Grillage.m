%% Anm_Grillage Class
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <anm.html *Anm*> in the <main.html LESM> (Linear Elements
% Structure Model) program. This sub-class implements abstract methods,
% declared in super-class *Anm*, that deal with grillage analysis model
% of linear structure elements.
%
%% Grillage Analysis
% A grillage model is a common form of analysis model for building
% stories and bridge decks. Its key features are:
%%%
% * It is a 2D model, which in LESM is considered in the XY-plane.
% * Beam elements are laid out in a grid pattern in a single plane,
%   rigidly connected at nodes. However, a grillage element might have
%   a hinge (rotation liberation) at an end or hinges at both ends.
% * It is assumed that a hinge in a grillage element releases continuity of
%   both bending and torsion rotations.
% * By assumption, there is only out-of-plane behavior, which
%   includes displacement transversal to the grillage plane, and
%   rotations about in-plane axes.
% * Internal forces at any cross-section of a grillage element are:
%   shear force (local z direction), bending moment (about local y direction),
%   and torsion moment (about local x direction).
%   By assumption, there is no axial force in a grillage element.
% * Each node of a grillage model has three d.o.f.'s: a transversal
%   displacement in Z direction, and rotations about the X and Y directions.
%
%% Class definition
% Definition of sub-class *Anm_Grillage* derived from super-class <anm.html Anm>.
classdef Anm_Grillage < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Grillage()
            include_constants;
            anm = anm@Anm(GRILLAGE_ANALYSIS,3);
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
                
                % Check for fixed rotation about global X direction
                if drv.nodes(n).ebc(4) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(1,n) = 1;
                end
                
                % Check for fixed rotation about global Y direction
                if drv.nodes(n).ebc(5) ~= FREE_DOF
                    drv.neqfixed = drv.neqfixed + 1;
                    drv.ID(2,n) = 1;
                end
                
                % Check for fixed translation in global Z direction
                if drv.nodes(n).ebc(3) ~= FREE_DOF
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
                    
                    % Add prescribed rotation about global X direction
                    id = drv.ID(1,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(4) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(4);
                    end
                    
                    % Add prescribed rotation about global Y direction
                    id = drv.ID(2,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(5) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(5);
                    end
                    
                    % Add prescribed displacement in global Z direction
                    id = drv.ID(3,n);
                    if (id > drv.neqfree) && (drv.nodes(n).prescDispl(3) ~= 0)
                        drv.D(id) = drv.nodes(n).prescDispl(3);
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
            % Compute torsion stiffness coefficients
            ket = elem.torsionStiffCoeff();
            
            % Compute flexural stiffness coefficients
            kef = elem.flexuralStiffCoeff_XZ();
            
            % Assemble element stiffness matrix in local system
            kel = [ ket(1,1)  0         0         ket(1,2)  0         0;
                    0         kef(2,2)  kef(2,1)  0         kef(2,4)  kef(2,3);
                    0         kef(1,2)  kef(1,1)  0         kef(1,4)  kef(1,3);
                    ket(2,1)  0         0         ket(2,2)  0         0;
                    0         kef(4,2)  kef(4,1)  0         kef(4,4)  kef(4,3);
                    0         kef(3,2)  kef(3,1)  0         kef(3,4)  kef(3,3) ];
        end
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  drv: handle to an object of the Drv class
        function nodalLoads(~,drv)
            for n = 1:drv.nnp
                if ~isempty(drv.nodes(n).nodalLoad)
                    
                    % Add applied moment about global X direction
                    id = drv.ID(1,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(4);
                    
                    % Add applied moment about global Y direction
                    id = drv.ID(2,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(5);
                    
                    % Add applied force in global Z direction
                    id = drv.ID(3,n);
                    drv.F(id) = drv.F(id) + drv.nodes(n).nodalLoad(3);
                    
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
            % Compute flexural (transversal) fixed end force components
            fef = load.flexuralDistribLoadFEF_XZ();
            
            % Assemble element fixed end force (FEF) vector in local system
            fel = [ 0;
                    fef(2);
                    fef(1);
                    0;
                    fef(4);
                    fef(3) ];
        end
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied thermal load (temperature variation).
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        function fel = elemLocThermalLoadFEF(~,load)
            % Compute flexural (transversal) fixed end force components
            fef = load.flexuralThermalLoadFEF_XZ();
            
            % Assemble element fixed end force (FEF) vector in local system
            fel = [ 0;
                    fef(2);
                    fef(1);
                    0;
                    fef(4);
                    fef(3) ];
        end
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        %  torsion_moment(1,2)
        %   Ti = torsion_moment(1) - init value
        %   Tf = torsion_moment(2) - final value
        %  bending_moment(1,2)
        %   Mi = bending_moment(1) - init value
        %   Mf = bending_moment(2) - final value
        %  shear_force(1,2)
        %   Qi = shear_force(1) - init value
        %   Qf = shear_force(2) - final value
        % Input arguments:
        %  elem: handle to an object of the Elem class
        function initIntForce(~,elem)
            elem.torsion_moment   = zeros(1,2);
            elem.bending_moment_Y = zeros(1,2);
            elem.shear_force_Z    = zeros(1,2);
        end
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        function assembleIntForce(~,elem,fel)
            elem.torsion_moment(1)   = elem.torsion_moment(1)   + fel(1);
            elem.torsion_moment(2)   = elem.torsion_moment(2)   + fel(4);
            elem.bending_moment_Y(1) = elem.bending_moment_Y(1) + fel(2);
            elem.bending_moment_Y(2) = elem.bending_moment_Y(2) + fel(5);
            elem.shear_force_Z(1)    = elem.shear_force_Z(1)    + fel(3);
            elem.shear_force_Z(2)    = elem.shear_force_Z(2)    + fel(6);
        end
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        %  intDispl(1,:) -> du (axial displacement)
        %  intDispl(2,:) -> dw (transversal displacement in local z-axis)
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
            % Compute transversal displacement shape functions vector
            Nw = elem.flexuralDisplShapeFcnVector_XZ(x);
            
            % Assemble displacement shape function matrix
            N = [ 0      0      0      0      0      0;
                  0      Nw(2)  Nw(1)  0      Nw(4)  Nw(3) ];
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements vector in local system,
        % in a given cross-section position, for the local analysis from
        % element loads (distributed loads and thermal loads).
        % Output:
        %  del: a 2x1 vector of with element internal displacements:
        %      du -> axial displacement (always zero)
        %      dw -> transversal displacement in local z-axis direction
        %      x  -> cross-section position on element local x-axis
        %  del = [ du(x);
        %          dw(x) ]
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        function del = lclAnlIntDispl(~,elem,x)
            % Initialize displacements vector resulting from local analysis
            del = [ 0;
                    0 ];
            
            % Add the contribution of transversal displacements resulting
            % from distributed loads (there is no axial displacement in a
            % grillage model)
            if (~isempty(elem.load.uniformLcl)) || (~isempty(elem.load.linearLcl))
                del(2) = elem.load.flexuralDistribLoadDispl_XZ(x);
            end
            
            % Add the contribution of transversal displacements resulting
            % from thermal loads (there is no axial displacement in a
            % grillage model)
            if (elem.load.tempVar_X ~= 0) || (elem.load.tempVar_Z ~= 0)
                del(2) = del(2) + elem.load.flexuralThermalLoadDispl_XZ(x);
            end
        end
    end
end