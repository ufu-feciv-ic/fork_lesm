%% Node Class
%
% This is a class in the Object Oriented Programming (OOP) paradigm
% that defines node objects in the <main.html LESM> (Linear Elements
% Structure Model) program.
%
% A node is a joint between two or more elements, or any element end, used
% to discretize the model. It is always considered as a tri-dimensional
% entity which may have applied loads or prescribed displacements towards
% one of its six degrees of freedom.
%
%% Class definition
% Definition of super-class *Node* as a <https://www.mathworks.com/help/matlab/ref/handle-class.html *handle*> class.
classdef Node < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id = 0;            % identification number
        coord = [];        % vector of coordinates on global system [X Y Z]
        ebc = [];          % vector of essential boundary condition flags [dx dy dz rx ry rz]
        nodalLoad = [];    % vector of applied load components [fx fy fz mx my mz]
        prescDispl = [];   % vector of prescribed displacement values [dx dy dz rx ry rz]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function node = Node(id,coord,ebc,nodalLoad,presDispl)
            if (nargin > 0)
                node.id = id;
                node.coord = coord;
                node.ebc = ebc;
                node.nodalLoad = nodalLoad;
                node.prescDispl = presDispl;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Counts total number of elements and number of hinged elements
        % connected to a node.
        % Output:
        %  tot: total number of elements connected to a node
        %  hng: number of hinged elements connected to a node
        % Input arguments:
        %  drv: handle to an object of the Drv class
        function [tot,hng] = elemsIncidence(node,drv)
            include_constants;
            
            % Initialize values
            tot = 0;
            hng = 0;
            
            for e = 1:drv.nel
                % Check if initial node of current element is the target one
                if drv.elems(e).nodes(1).id == node.id
                    tot = tot + 1;
                    if drv.elems(e).hingei == HINGED_END
                        hng = hng + 1;
                    end
                
                % Check if final node of current element is the target one
                elseif drv.elems(e).nodes(2).id == node.id
                    tot = tot + 1;
                    if drv.elems(e).hingef == HINGED_END
                        hng = hng + 1;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a Node object.
        function clean(node)
            node.id = 0;
            node.coord = [];
            node.ebc = [];
            node.nodalLoad = [];
            node.prescDispl = [];
        end
    end
end