%% Lelem (Load Element) Class
%
% This is an abstract super-class in the Object Oriented Programming (OOP)
% paradigm that generically specifies an element load object in the <main.html
% LESM> (Linear Elements Structure Model) program.
%
% This super-class is abstract because one cannot instantiate objects from
% it, since it contains abstract methods that should be implemented in
% (derived) sub-classes.
%
% Essentially, this super-class declares abstract methods that define
% the generic flexural behavior of linear elements subjected to internal loads
% in the <main.html LESM> (Linear Elements Structure Model) program.
% These abstract methods are the functions that should be implemented in a
% (derived) sub-class that deals with specific flexural behavior, such as Navier
% (Euler-Bernoulli) or Timoshenko flexural behavior.
% In addition, this super-class also implements methods that consider
% everything that does not depend on flexural behavior, i.e., methods that
% deal with axial load.
%
% An element load object is responsible for defining the response of a
% linear element to loading effects, such as the computation of element fixed
% end forces (FEF), element equivalent nodal loads (ENL), and element
% internal displacements. There is a mutual relation between an object of
% the *Lelem* class and an object of the <elem.html *Elem*> class, since each
% element has its own load properties and each load object is associated
% with one element.
%
%% Element Load Types
%
% There are three types of element loads considered in the LESM program:
%
% * Uniformely distributed force on elements, spanning its entire
%   length, in local or in global axes directions.
% * Linearly distributed force on elements, spanning its entire length,
%   in local or in global axes directions.
% * Temperature variation along elements local axes (temperature gradient).
%
% It is assumed that there is a single load case.
%
%% Sub-classes
% To handle different types of flexural behavior, there are two sub-classes
% derived from this super-class:
%%%
% * <lelem_navier.html Navier (Euler-Bernoulli) flexural behavior>.
% * <lelem_timoshenko.html Timoshenko flexural behavior>.
%
% In truss models, the two types of elements may be used indistinguishably,
% since there is no bending behavior of a truss element, and
% Euler-Bernoulli elements and Timoshenko elements are equivalent for the
% axial behavior.
%
%% Class definition
% Definition of super-class *Lelem* as a <https://www.mathworks.com/help/matlab/ref/handle-class.html *handle*> class.
classdef Lelem < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        elem = [];           % handle to an object of the Elem class
        uniformDir = 0;      % flag for uniform load direction
        uniformGbl = [];     % vector of uniformly distributed load components in global system [qx qy qz]
        uniformLcl = [];     % vector of uniformly distributed load components in local system  [qx qy qz]
        linearDir = 0;       % flag for linear load direction
        linearGbl = [];      % vector of linearly distributed load components in global system [qxi qyi qzi qxf qyf qzf]
        linearLcl = [];      % vector of linearly distributed load components in local system  [qxi qyi qzi qxf qyf qzf]
        tempVar_X = 0;       % temperature variation on element center of gravity
        tempVar_Y = 0;       % temperature gradient relative to local y-axis
        tempVar_Z = 0;       % temperature gradient relative to local z-axis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function load = Lelem(elem)
            load.elem = elem;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xy-plane for an applied linearly distributed load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        fef = flexuralDistribLoadFEF_XY(load)
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xz-plane for an applied linearly distributed load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        fef = flexuralDistribLoadFEF_XZ(load)
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xy-plane for an applied thermal load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        fef = flexuralThermalLoadFEF_XY(load)
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xz-plane for an applied thermal load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        fef = flexuralThermalLoadFEF_XZ(load)
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xy-plane at
        % given cross-section position for an applied linearly distributed load.
        % Output:
        %  dv: transversal displacement in local y-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        dv = flexuralDistribLoadDispl_XY(load,x)
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xz-plane at
        % given cross-section position for an applied linearly distributed load.
        % Output:
        %  dw: transversal displacement in local z-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        dw = flexuralDistribLoadDispl_XZ(load,x)
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xy-plane at
        % given cross-section position for an applied thermal load.
        % Output:
        %  dv: transversal displacement in local y-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        dv = flexuralThermalLoadDispl_XY(load,x)
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xz-plane at
        % given cross-section position for an applied thermal load.
        % Output:
        %  dw: transversal displacement in local z-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        dw = flexuralThermalLoadDispl_XZ(load,x)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Sets uniformly distributed load in global and local system.
        % Input arguments:
        %  unifload: vector of uniformly distributed load components
        %  dir:      specified direction of the distributed load components
        function setUnifLoad(load,unifload,dir)
            include_constants
            rot = load.elem.T;
            
            if dir == GLOBAL_LOAD
                load.uniformGbl = unifload;
                load.uniformLcl = rot * unifload';
                
            elseif dir == LOCAL_LOAD
                load.uniformLcl = unifload;
                load.uniformGbl = rot' * unifload';
            end
        end
        
        %------------------------------------------------------------------
        % Sets linearly distributed load in global and local system.
        % Input arguments:
        %  linearload: vector of linearly distributed load components
        %  dir:        specified direction of the distributed load components
        function setLinearLoad(load,linearload,dir)
            include_constants
            rot = blkdiag(load.elem.T,load.elem.T);
            
            if dir == GLOBAL_LOAD
                load.linearGbl = linearload;
                load.linearLcl = rot * linearload';
                
            elseif dir == LOCAL_LOAD
                load.linearLcl = linearload;
                load.linearGbl = rot' * linearload';
            end
        end
        
        %------------------------------------------------------------------
        % Generates element axial fixed end force (FEF) vector for an
        % applied linearly distributed load.
        % Output:
        %  fea: 2 position vector with axial FEF's:
        %       fea(1) -> axial force at initial node
        %       fea(2) -> axial force at final node
        function fea = axialDistribLoadFEF(load)
            % Initialize axial load values at end nodes
            qxi = 0;
            qxf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qxi = qxi + load.uniformLcl(1);
                qxf = qxi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qxi = qxi + load.linearLcl(1);
                qxf = qxf + load.linearLcl(4);
            end
            
            % Check if axial load is not null over element
            if (qxi ~= 0) || (qxf ~= 0)
                L = load.elem.length;
                
                % Separate uniform portion from linear portion of axial load
                qx0 = qxi;
                qx1 = qxf - qxi;
                
                % Calculate fixed end force vector
                fea = [ -(qx0*L/2 + qx1*L/6);
                        -(qx0*L/2 + qx1*L/3) ];
               
            else
                
                fea = [ 0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element axial fixed end force (FEF) vector for an
        % applied thermal load.
        % Output:
        %  fea: 2 position vector with axial FEF's:
        %       fea(1,1) -> axial force at initial node
        %       fea(2,1) -> axial force at final node
        function fea = axialThermalLoadFEF(load)
            % Get temperature variation on element center of gravity
            dtx = load.tempVar_X;
            
            % Check if temperature variation is not null
            if dtx ~= 0
                E     = load.elem.material.elasticity;
                alpha = load.elem.material.thermExp;
                A     = load.elem.section.area_x;
                
                % Calculate fixed end force vector
                fea = [ E * A * alpha * dtx;
                       -E * A * alpha * dtx ];
            
            else
                
                fea = [ 0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Computes element axial displacement at given cross-section position
        % for an applied linearly distributed load assuming a fixed end
        % condition.
        % Output:
        %  du: axial displacement in longitudinal direction
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function du = axialDistribLoadDispl(load,x)
            include_constants;
            
            % Initialize axial load values at end nodes
            qxi = 0;
            qxf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qxi = qxi + load.uniformLcl(1);
                qxf = qxi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qxi = qxi + load.linearLcl(1);
                qxf = qxf + load.linearLcl(4);
            end
            
            % Check if transversal load is not null over element
            if (qxi ~= 0) || (qxf ~= 0)
                E = load.elem.material.elasticity;
                A = load.elem.section.area_x;
                L = load.elem.length;

                x2 = x^2;
                x3 = x2 * x;
                EA = E  * A;

                % Separate uniform portion from linear portion of axial load
                qx0 = qxi;
                qx1 = qxf - qxi;

                % Calculate axial displacements from uniform load portion (up0)
                % and from linear axial load portion (up1)
                up0 = qx0/EA * (L*x/2 - x2/2);
                up1 = qx1/EA * (L*x/6 - x3/(6*L));
                
                du = up0 + up1;
                
            else
                du = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element axial displacement at given cross-section position
        % for an applied thermal load assuming a fixed end condition.
        % Output:
        %  du: axial displacement in longitudinal direction
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function du = axialThermalLoadDispl(~,~)
            % Axial displacement is null when an element with fixed ends is
            % subjected to a temperature variation
            du = 0;
        end
        
        %------------------------------------------------------------------
        % Computes element equivalent nodal load (ENL) vector in global
        % system for an applied distributed load.
        % Output:
        %  feg: equivalent nodal load vector in global system
        function feg = gblDistribLoadENL(load)
            % Compute and store element fixed end forces (FEF) in local
            % system for a distributed load
            load.elem.fel_distribLoad = load.elem.anm.elemLocDistribLoadFEF(load);
            
            % Transform element fixed end forces in local system to
            % equivalent nodal loads in global system
            feg = -load.elem.rot' * load.elem.fel_distribLoad;
        end
        
        %------------------------------------------------------------------
        % Computes element equivalent nodal load (ENL) vector in global
        % system for an applied thermal load.
        % Output:
        %  feg: equivalent nodal load vector in global system
        function feg = gblThermalLoadENL(load)
            % Compute and store element fixed end forces (FEF) in local
            % system for a thermal load
            load.elem.fel_thermalLoad = load.elem.anm.elemLocThermalLoadFEF(load);
            
            % Transform element fixed end forces in local system to
            % equivalent nodal loads in global system
            feg = -load.elem.rot' * load.elem.fel_thermalLoad;
        end
        
        %------------------------------------------------------------------
        % Clean data structure of a Lelem object.
        function clean(load)
            load.elem = [];
            load.uniformDir = 0;
            load.uniformGbl = [];
            load.uniformLcl = [];
            load.linearDir = 0;
            load.linearGbl = [];
            load.linearLcl = [];
            load.tempVar_X = 0;
            load.tempVar_Y = 0;
            load.tempVar_Z = 0;
        end
    end
end