%% Anm (Analysis Model) Class
%
% This is an abstract super-class in the Object Oriented Programming (OOP)
% paradigm that generically specifies an analysis model in the <main.html
% LESM> (Linear Elements Structure Model) program.
%
% This super-class is abstract because one cannot instantiate objects from
% it, since it contains abstract methods that should be implemented in
% (derived) sub-classes.
%
% Essentially, this super-class declares abstract methods that define
% the generic behavior of an analysis model in the <main.html
% LESM> (Linear Elements Structure Model) program. These abstract methods
% are the functions that should be implemented in a (derived) sub-class
% that deals with specific types of analysis, such as a 2D truss, 2D frame,
% grillage, 3D truss or a 3D frame model.
%
% An object of the *Anm* class "projects" the generic 3D behavior of linear
% element objects of the <elem.html *Elem*> class and node objects of the
% <node.html *Node*> class to a specific model behavior.
% The implementations of an abstract method in each sub-class differ in the
% use of the degrees of freedom of each analysis model.
% Only one analysis model object is created during the analysis process.
%
%% Sub-classes
% To handle different types of analysis models, there are five sub-classes
% derived from this super-class:
%%%
% * <anm_truss2d.html 2D truss analysis model>.
% * <anm_frame2d.html 2D frame analysis model>.
% * <anm_grillage.html Grillage analysis model>.
% * <anm_truss3d.html 3D truss analysis model>.
% * <anm_frame3d.html 3D frame analysis model>.
%
%% Class definition
% Definition of super-class *Anm* as a <https://www.mathworks.com/help/matlab/ref/handle-class.html *handle*> class.
classdef Anm < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        analysis_type = 0;   % flag for type of analysis model
        ndof = 0;            % number of degrees-of-freedom per node
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm(type,ndof)
            anm.analysis_type = type;
            anm.ndof = ndof;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Assembles element d.o.f. (degree of freedom) rotation transformation
        % matrix from global system to local system.
        % Output:
        %  rot: rotation transformation matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        rot = gblToLocElemRotMtx(~,elem)
        
        %------------------------------------------------------------------
        % Initializes global d.o.f (degree of freedom) numbering ID matrix
        % with ones and zeros, and counts total number of equations of free
        % d.o.f.'s and total number of equations of fixed d.o.f.'s.
        %  ID matrix initialization:
        %  if ID(k,n) = 0, d.o.f. k of node n is free.
        %  if ID(k,n) = 1, d.o.f. k of node n is constrained by support.
        % Input arguments:
        %  drv: handle to an object of the Drv class
        setupDOFNum(anm,drv)
        
        %------------------------------------------------------------------
        % Adds prescribed displacements (known support settlement values)
        % to global displacement vector.
        % Avoids storing a prescribed displacement component in a position
        % of the global displacement vector that corresponds to a free d.o.f.
        % Input arguments:
        %  drv: handle to an object of the Drv class
        setupPrescDispl(~,drv)
        
        %------------------------------------------------------------------
        % Assembles element stiffness matrix in local system.
        % Output:
        %  kel: target element stiffness matrix in local system
        % Input arguments:
        %  elem: handle to an object of the Elem class
        kel = elemLocStiffMtx(~,elem)
        
        %------------------------------------------------------------------
        % Adds nodal load components to global forcing vector,
        % including the terms that correspond to constrained d.o.f.
        % Input arguments:
        %  drv: handle to an object of the Drv class
        nodalLoads(~,drv)
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied distributed load.
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        fel = elemLocDistribLoadFEF(~,load)
        
        %------------------------------------------------------------------
        % Assembles element fixed end force (FEF) vector in local system
        % for an applied thermal load (temperature variation).
        % Output:
        %  fel: element fixed end force vector in local system
        % Input arguments:
        %  load: handle to an object of the Lelem class
        fel = elemLocThermalLoadFEF(~,load)
        
        %------------------------------------------------------------------
        % Initializes element internal forces arrays with null values.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        initIntForce(~,elem)
        
        %------------------------------------------------------------------
        % Assembles contribution of a given internal force vector to
        % element arrays of internal forces.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  fel: element internal force vector in local system
        assembleIntForce(~,elem,fel)
        
        %------------------------------------------------------------------
        % Initializes element internal displacements array with null values.
        % Each element is discretized in 50 cross-sections, where internal
        % displacements are computed.
        % Input arguments:
        %  elem: handle to an object of the Elem class
        initIntDispl(~,elem)
        
        %------------------------------------------------------------------
        % Assembles displacement shape functions matrix evaluated at a 
        % given cross-section position.
        % Output:
        %  N: displacement shape function matrix
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position
        N = displShapeFcnMtx(~,elem,x)
        
        %------------------------------------------------------------------
        % Computes element internal displacements vector in local system,
        % in a given cross-section position, for the local analysis from
        % element loads (distributed loads and thermal loads).
        % Output:
        %  del: element internal displacements vector in a given
        %       cross-section position
        % Input arguments:
        %  elem: handle to an object of the Elem class
        %  x: cross-section position on element local x-axis
        del = lclAnlIntDispl(anm,elem,x)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of an Anm object.
        function clean(anm)
            anm.analysis_type = 0;
            anm.ndof = 0;
        end
    end
end