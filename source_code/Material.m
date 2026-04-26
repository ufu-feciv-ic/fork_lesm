%% Material Class
%
% This is a class in the Object Oriented Programming (OOP) paradigm
% that defines material objects in the <main.html LESM> (Linear Elements
% Structure Model) program.
%
% All materials in LESM are considered to have linear elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in
% all directions.
%
%% Class definition
% Definition of super-class *Material* as a <https://www.mathworks.com/help/matlab/ref/handle-class.html *handle*> class.
classdef Material < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 0;           % identification number
        elasticity = 0;   % elasticity modulus
        poisson = 0;      % poisson ratio
        shear = 0;        % shear modulus
        thermExp = 0;     % thermal expansion coefficient
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function material = Material(id,e,v,te)
            if (nargin > 0)
                material.id = id;
                material.elasticity = e;
                material.poisson = v;
                material.shear = e / (2 * (1 + v));
                material.thermExp = te;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Cleans data structure of a Material object.
        function clean(material)
            material.id = 0;
            material.elasticity = 0;
            material.poisson = 0;
            material.shear = 0;
            material.thermExp = 0;
        end
    end
end