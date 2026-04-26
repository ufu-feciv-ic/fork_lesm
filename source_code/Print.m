%% Print Class
%
% This is an abstract super-class in the Object Oriented Programming (OOP)
% paradigm that generically specifies a printing object in the <main.html LESM>
% (Linear Elements Structure Model) program.
%
% This super-class is abstract because one cannot instantiate objects from
% it, since it contains abstract methods that should be implemented in
% (derived) sub-classes.
%
% Essentially, this super-class declares abstract methods and define
% public methods that print model information and analysis results of the
% <main.html LESM> (Linear Elements Structure Model) program in a text
% file or in the default output (MATLAB command window).
% These abstract methods are the functions that should be implemented
% in a (derived) sub-class that deals with specific types of models,
% such as a 2D truss, 2D frame, grillage, 3D truss or a 3D frame model.
%
%% Sub-classes
% To handle different types of analysis models, there are five sub-classes
% derived from this super-class:
%%%
% * <print_truss2d.html Print 2D truss model and results>.
% * <print_frame2d.html Print 2D frame model and results>.
% * <print_grillage.html Print grillage model and results>.
% * <print_truss3d.html Print 3D truss model and results>.
% * <print_frame3d.html Print 3D frame model and results>.
%
%% Class definition
% Definition of super-class *Print* as a <https://www.mathworks.com/help/matlab/ref/handle-class.html *handle*> class.
classdef Print < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        txt = 0;    % output identifier (1 = MATLAB command window)
        drv = [];   % handle to an object of the Drv class
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function print = Print(txt)
            print.txt = txt;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Prints analyis results.
        results(print)
        
        %------------------------------------------------------------------
        % Prints analyis model type.
        analysisLabel(print)
        
        %------------------------------------------------------------------
        % Prints model description.
        % Output:
        %  n_nl: number of nodes with applied loads
        %  n_pd: number of nodes with prescribed displacement
        %  n_ul: number of elements with uniformly distributed loads
        %  n_ll: number of elements with linearly distributed loads
        %  n_tv: number of elements with temperature variation
        [n_nl,n_pd,n_ul,n_ll,n_tv] = modelDescrip(print)
        
        %------------------------------------------------------------------
        % Prints cross-section properties.
        section(print)
        
        %------------------------------------------------------------------
        % Prints nodal support conditions.
        nodalSupport(print)
        
        %------------------------------------------------------------------
        % Prints nodal loads.
        % Input arguments:
        %  n_nl: number of nodes with applied loads
        nodalLoads(print,n_nl)
        
        %------------------------------------------------------------------
        % Prints nodal prescribed displacements.
        % Input arguments:
        %  n_pd: number of nodes with prescribed displacement
        nodalPrescDisp(print,n_pd)
        
        %------------------------------------------------------------------
        % Prints elements information.
        elements(print)
        
        %------------------------------------------------------------------
        % Prints uniformly distributed loads information.
        % Input arguments:
        %  n_ul: number of elements with uniformly distributed loads
        unifElementLoads(print,n_ul)
        
        %------------------------------------------------------------------
        % Prints linearly distributed loads information.
        % Input arguments:
        %  n_ll: number of elements with linearly distributed loads
        linearElementLoads(print,n_ll)
        
        %------------------------------------------------------------------
        % Prints thermal loads information.
        % Input arguments:
        %  n_tv: number of elements with temperature variation
        temperatureVariation(print,n_tv)
        
        %------------------------------------------------------------------
        % Prints results of nodal displacement/rotation.
        nodalDisplRot(print)
        
        %------------------------------------------------------------------
        % Prints results of support reactions.
        reactions(print)
        
        %------------------------------------------------------------------
        % Prints results of internal forces at element nodes.
        intForces(print)
        
        %------------------------------------------------------------------
        % Prints results of elements internal displacements.
        elemDispl(print)
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Prints header of analysis results.
        function header(print)
            fprintf(print.txt, '\n=========================================================\n');
            fprintf(print.txt, ' LESM - Linear Elements Structure Model analysis program\n');
            fprintf(print.txt, '    PONTIFICAL CATHOLIC UNIVERSITY OF RIO DE JANEIRO\n');
            fprintf(print.txt, '    DEPARTMENT OF CIVIL AND ENVIRONMENTAL ENGINEERING\n');
            fprintf(print.txt, '                          AND\n');
            fprintf(print.txt, '               TECGRAF/PUC-RIO INSTITUTE\n');
            fprintf(print.txt, '=========================================================\n');
        end
        
        %------------------------------------------------------------------
        % Prints material properties.
        function material(print)
            fprintf(print.txt, '\n\n-------------------------------------\n');
            fprintf(print.txt, 'M A T E R I A L  P R O P E R T I E S\n');
            fprintf(print.txt, '-------------------------------------\n');
            if print.drv.nmat > 0
                fprintf(print.txt, ' MATERIAL     E [MPa]      G [MPa]     POISSON      THERMAL EXP. COEFF. [/°C]\n');
                
                for m = 1:print.drv.nmat
                    fprintf(print.txt, '%4d   %13.0f   %10.0f   %8.2f   %20.3d\n', m,...
                            1e-3*print.drv.materials(m).elasticity,...
                            1e-3*print.drv.materials(m).shear,...
                            print.drv.materials(m).poisson,...
                            print.drv.materials(m).thermExp);
                end
            else
                fprintf(print.txt, ' NO MATERIAL\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal coordinates.
        function nodalCoords(print)
            fprintf(print.txt, '\n\n--------------------------------\n');
            fprintf(print.txt, 'N O D A L  C O O R D I N A T E S \n');
            fprintf(print.txt, '--------------------------------\n');
            fprintf(print.txt, ' NODE     COORD X [m]      COORD Y [m]      COORD Z [m]    \n');
            
            for n = 1:print.drv.nnp
                fprintf(print.txt, '%4d   %11.3f   %14.3f   %14.3f\n', n,...
                        print.drv.nodes(n).coord(1),...
                        print.drv.nodes(n).coord(2),...
                        print.drv.nodes(n).coord(3));
            end
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a Print object.
        function clean(print)
            print.txt = [];
            print.drv = [];
        end
    end
end