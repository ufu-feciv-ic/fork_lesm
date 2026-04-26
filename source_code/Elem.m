%% Elem (Element) Class
%
% This is an abstract super-class in the Object Oriented Programming (OOP)
% paradigm that generically specifies a linear element in the <main.html
% LESM> (Linear Elements Structure Model) program.
%
% This super-class is abstract because one cannot instantiate objects from
% it, since it contains abstract methods that should be implemented in
% (derived) sub-classes.
%
% Essentially, this super-class declares abstract methods that define
% the generic flexural behavior of a linear element in the <main.html
% LESM> (Linear Elements Structure Model) program. These abstract methods
% are the functions that should be implemented in a (derived) sub-class
% that deals with specific flexural behavior, such as Navier (Euler-Bernoulli)
% or Timoshenko flexural behavior.
% In addition, this super-class also implements methods that consider
% everything that does not depend on flexural behavior, i.e., methods that
% deal with axial and torsional behavior.
%
% The *Elem* class generically handles a three-dimensional behavior of a
% linear element. An object of the <anm.html *Anm*> class is responsible
% to "project" this generic 3D behavior to a specific model behavior,
% such as a 2D frame, 2D truss, grillage, 3D truss or a 3D frame model.
% Therefore, one of the properties of an object of the *Elem* class is a
% handle to a target <anm.html *Anm*> object.
%
%% Element Local Coordinate System
%
% In 2D models of LESM, the local axes of an element are defined
% uniquely in the following manner:
%
% * The local z-axis of an element is always in the direction of the
%   global Z-axis, which is perpendicular to the model plane and its
%   positive direction points out of the screen.
% * The local x-axis of an element is its longitudinal axis, from its
%   initial node to its final node.
% * The local y-axis of an element lays on the global XY-plane and is
%   perpendicular to the element x-axis in such a way that the
%   cross-product (x-axis * y-axis) results in a vector in the
%   global Z direction.
%
% In 3D models of LESM, the local y-axis and z-axis are defined by an
% auxiliary vector vz = (vzx,vzy,vzz), which is an element property and
% should be specified as an input data of each element:
%
% * The local x-axis of an element is its longitudinal axis, from its
%   initial node to its final node.
% * The auxiliary vector lays in the local xz-plane of an element, and the
%   cross-product (vz * x-axis) defines the the local y-axis vector.
% * The direction of the local z-axis is then calculated with the
%   cross-product (x-axis * y-axis).
%
% In 2D models, the auxiliary vector vz is automatically set to (0,0,1).
% In 3D models, it is important that the auxiliary vector is not parallel
% to the local x-axis; otherwise, the cross-product (vz * x-axis) is zero
% and local y-axis and z-axis will not be defined.
%
%% Sub-classes
% To handle different types of flexural behavior, there are two sub-classes
% derived from this super-class:
%%%
% * <elem_navier.html Navier (Euler-Bernoulli) flexural behavior>.
% * <elem_timoshenko.html Timoshenko flexural behavior>.
%
% In truss models, the two types of elements may be used indistinguishably,
% since there is no bending behavior of a truss element, and
% Euler-Bernoulli elements and Timoshenko elements are equivalent for the
% axial behavior.
%
%% Class definition
% Definition of super-class *Elem* as a <https://www.mathworks.com/help/matlab/ref/handle-class.html *handle*> class.
classdef Elem < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        nen = 2;                 % number of nodes (always equal to 2)
        type = 0;                % flag for type of element
        anm = [];                % handle to an object of the Anm class
        load = [];               % handle to an object of the Lelem class
        material = [];           % handle to an object of the Material class
        section = [];            % handle to an object of the Section class
        nodes = [];              % vector of handles to objects of the Node class [initial_node final_node]
        hingei = 0;              % flag for rotation liberation at initial node
        hingef = 0;              % flag for rotation liberation at final node
        length = 0;              % element length
        cosine_X = 0;            % cosine of orientation angle with global X-axis
        cosine_Y = 0;            % cosine of orientation angle with global Y-axis
        cosine_Z = 0;            % cosine of orientation angle with global Z-axis
        vz = [];                 % auxiliary vector in local xz-plane [vz_X vz_Y vz_Z]
        gle = [];                % gather vector (stores element d.o.f. eqn. numbers)
        T = [];                  % basis rotation transformation matrix between two coordinate systems
        rot = [];                % rotation transformation matrix from global system to local system
        kel = [];                % stiffness matrix in local system
        fel_distribLoad = [];    % fixed end forces in local system for a distributed load
        fel_thermalLoad = [];    % fixed end forces in local system for a thermal load
        axial_force = [];        % array of axial internal forces at element ends
        shear_force_Y = [];      % array of shear internal forces at element ends in local y-axis
        shear_force_Z = [];      % array of shear internal forces at element ends in local z-axis
        bending_moment_Y = [];   % array of bending internal moments at element ends in local y-axis
        bending_moment_Z = [];   % array of bending internal moments at element ends in local z-axis
        torsion_moment = [];     % array of torsional internal moments at element ends
        intDispl = [];           % array of axial and transversal internal displacements
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Elem(type)
            elem.type = type;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
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
        kef = flexuralStiffCoeff_XY(elem)
        
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
        kef = flexuralStiffCoeff_XZ(elem)
        
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
        Nv = flexuralDisplShapeFcnVector_XY(elem,x)
        
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
        Nw = flexuralDisplShapeFcnVector_XZ(elem,x)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes basis rotation transformation matrix between local and
        % global coordinate systems of an element.
        % Output:
        %  T: basis rotation transformation matrix
        function T = rotTransMtx(elem)
            % Get nodal coordinates
            xi = elem.nodes(1).coord(1);
            yi = elem.nodes(1).coord(2);
            zi = elem.nodes(1).coord(3);
            xf = elem.nodes(2).coord(1);
            yf = elem.nodes(2).coord(2);
            zf = elem.nodes(2).coord(3);
            
            % Calculate element local x-axis
            x = [xf-xi, yf-yi, zf-zi];
            x = x / norm(x);
            
            % Calculate element local y-axis
            y = cross(elem.vz,x);
            y = y / norm(y);
            
            % Calculate element local z-axis
            z = cross(x,y);
            
            % Global axes
            X = [1,0,0];
            Y = [0,1,0];
            Z = [0,0,1];
            
            % Calculate angle cosines between local x-axis and global axes
            cxX = dot(x,X);
            cxY = dot(x,Y);
            cxZ = dot(x,Z);
            
            % Calculate angle cosines between local y-axis and global axes
            cyX = dot(y,X);
            cyY = dot(y,Y);
            cyZ = dot(y,Z);
            
            % Calculate angle cosines between local z-axis and global axes
            czX = dot(z,X);
            czY = dot(z,Y);
            czZ = dot(z,Z);
            
            % Assemble basis transformation matrix
            T = [ cxX cxY cxZ;
                  cyX cyY cyZ;
                  czX czY czZ ];
        end
        
        %------------------------------------------------------------------
        % Computes element stiffness matrix in global system.
        % Output:
        %  keg: element stiffness matrix in global system
        function keg = gblStiffMtx(elem)
            % Compute and store element stiffness matrix in local system
            elem.kel = elem.anm.elemLocStiffMtx(elem);
            
            % Transform element stiffness matrix from local to global system
            keg = elem.rot' * elem.kel * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Generates element axial stiffness coefficient matrix.
        % Output:
        % kea: a 2x2 matrix with axial stiffness coefficients
        function kea = axialStiffCoeff(elem)
            E = elem.material.elasticity;
            A = elem.section.area_x;
            L = elem.length;
            
            kea = [ E*A/L  -E*A/L;
                   -E*A/L   E*A/L ];
        end
        
        %------------------------------------------------------------------
        % Generates element torsion stiffness coefficient matrix.
        % Output:
        %  ket: a 2x2 matrix with torsion stiffness coefficients
        function ket = torsionStiffCoeff(elem)
            include_constants;
            
            G  = elem.material.shear;
            Jt = elem.section.inertia_x;
            L  = elem.length;
            
            if (elem.hingei == CONTINUOUS_END) && (elem.hingef == CONTINUOUS_END)
                
                ket = [ G*Jt/L  -G*Jt/L;
                       -G*Jt/L   G*Jt/L ];
                
            else
                
                ket = [ 0  0;
                        0  0 ];
                
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal force vector in local system at end nodes,
        % for the global analysis (from nodal displacements and rotations).
        % Output:
        %  fel: element internal force vector in local system at end nodes
        % Input arguments:
        %  drv: handle to an object of the Drv class
        function fel = gblAnlIntForce(elem,drv)
            % Get nodal displacements and rotations at element end nodes
            % in global system
            dg = drv.D(elem.gle);
            
            % Transform solution from global system to local system
            dl = elem.rot * dg;
            
            % Compute internal force vector in local system
            fel = elem.kel * dl;
        end
        
        %------------------------------------------------------------------
        % Evaluates element axial displacement shape functions at a given
        % cross-section position of given element.
        % Axial shape functions give the longitudinal displacement of an
        % element subjected to unit nodal displacements in axial directions.
        % Output:
        %  Nu: a 1x2 vector with axial displacement shape functions:
        %      D -> axial displacement in local x-axis
        %      ini -> initial node
        %      end -> end node
        %   Nu = [ N_Dini  N_Dend ]
        % Input arguments:
        %  x: cross-section position on element local x-axis
        function Nu = axialDisplShapeFcnVector(elem,x)
            L = elem.length;
            
            Nu1 = 1 - x/L;
            Nu2 = x/L;
            
            Nu = [ Nu1  Nu2 ];
        end
        
        %------------------------------------------------------------------
        % Computes element axial and transversal internal displacements
        % vector in local system at a given cross-section position,
        % for the global analysis (from nodal displacements and rotations).
        % Output:
        %  del: element internal displacements vector at given
        %       cross-section position
        % Input arguments:
        %  drv: handle to an object of the Drv class
        %  x: cross-section position on element local x-axis
        function del = gblAnlIntDispl(elem,drv,x)
            % Get nodal displacements and rotations at element end nodes
            % in global system
            dg = drv.D(elem.gle);
            
            % Transform solution from global system to local system
            dl = elem.rot * dg;
            
            % Evaluate element shape function matrix at given cross-section
            N = elem.anm.displShapeFcnMtx(elem,x);
            
            % Compute internal displacements vector from global analysis
            del = N * dl;
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of an Elem object.
        function clean(elem)
            elem.nen = 2;
            elem.type = 0;
            elem.anm = [];
            elem.material = [];
            elem.section = [];
            elem.nodes = [];
            elem.hingei = 0;
            elem.hingef = 0;
            elem.length = 0;
            elem.cosine_X = 0;
            elem.cosine_Y = 0;
            elem.cosine_Z = 0;
            elem.vz = [];
            elem.gle = [];
            elem.rot = [];
            elem.T = [];
            elem.kel = [];
            elem.fel_distribLoad = [];
            elem.fel_thermalLoad = [];
            elem.load = [];
            elem.axial_force = [];
            elem.shear_force_Y = [];
            elem.shear_force_Z = [];
            elem.bending_moment_Y = [];
            elem.bending_moment_Z = [];
            elem.torsion_moment = [];
            elem.intDispl = [];
        end
    end
end