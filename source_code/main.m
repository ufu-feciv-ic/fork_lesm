%% LESM - Linear Elements Structure Model
%
% This is the main driver file of LESM. This is a MATLAB program for
% linear-elastic, displacement-based, static analysis of two-dimensional
% and three-dimensional linear elements structure models, using the direct
% stiffness method. For each structural analysis, the program assembles a
% system of equilibrium equations, solves the system and displays the
% analysis results.
%
% The program handles five different types of analysis models of reticulated
% structures: 2D (plane) Truss, 2D (plane) Frame, Grillage,
% 3D (spatial) Truss, and 3D (spatial) Frame.
%
% In addition, two types of flexural behavior for linear elements are
% considered: Euler-Bernoulli (Navier theory) elements and Timoshenko elements.
%
% The Object-Oriented Programming (OOP) paradigm is adopted in the
% implementation of the analysis process. The use of OOP is justified by
% the clarity, organization and generality that this programming paradigm
% provides to the code.
%
% The program may be used in a non-graphical version or in a GUI
% (Graphical User Interface) version.
% The non-graphical version reads a structural model from a neutral
% format file and prints model information and analysis results in the
% the default output (MATLAB command window).
% In the GUI version, a user may create a structural model with
% attributes through the program graphical interface. The program can
% save and read a structural model data stored in a neutral format file.
% The standalone executable of the GUI version can be downloaded from the
% LESM website: <https://web.tecgraf.puc-rio.br/lesm/>
%
%% Authors
%%%
% * Luiz Fernando Martha (lfm@tecgraf.puc-rio.br)
%%%
% * Rafael Lopez Rangel (rafaelrangel@tecgraf.puc-rio.br)
%%%
% Pontifical Catholic University of Rio de Janeiro - PUC-Rio
%
% Department of Civil and Environmental Engineering and Tecgraf Institute
% of Technical-Scientific Software Development of PUC-Rio (Tecgraf/PUC-Rio)
%
%% History
% @version 1.1
%
% Initial version: August 2017
%%%
% Initially prepared for the course ENG 1240 - Analise de Estruturas III,
% 2016, second term, Department of Civil Engineering, PUC-Rio.
% Developed in the undergraduate thesis "Development of a Graphic Program
% for Structural Analysis of Linear Element Models", by Rafael Lopez Rangel,
% advised by Professor Luiz Fernando Martha, Department of Civil Engineering,
% PUC-Rio, December, 2016.
% Extended in the undergraduate thesis "Extensao de Programa Grafico para
% Analise de Trelicas e Porticos Espaciais via MATLAB e GUI", by Murilo
% Felix Filho, advised by Professor Luiz Fernando Martha, Department of
% Civil Engineering, PUC-Rio, July, 2017.
%
% This version of the program is an accompanying software for the book
% "Analise Matricial de Estruturas com Orientacao a Objetos", Elsevier, 2018,
% by Luiz Fernando Martha.
%
%% Linear Elements Models
%
% Linear element models can be trusses, frames or grillages, which are made
% up of uniaxial members with one dimension much larger than the others,
% such as bars (members with axial behavior) or beams (members with flexural
% behavior). In the present context, members of any type are generically
% referred to as elements.
%
% The analysis process of this type of model is based on a nodal
% dicretization, where a node is a joint between two or more elements,
% or any element end.
%
% The final analysis result is composed by a the results from a global
% analysis and a local analysis. The result from the global analysis is
% computed with nodal displacements and rotations, and the result of a
% local analysis is computed with the effect of distributed loads over
% elements.
% These analyzes use analytical expressions of linear element behavior.
% Therefore, the sum of the two analysis provides the analytical solution
% of the problem.
%
%% Coordinate Systems
%
% Two types of coordinate systems are used to locate points or define
% directions in space. Both of these systems are cartesian, orthogonal and
% right-handed.
%%%
% *Global system*
%
% The global system XYZ is an absolute reference. Its origin is fixed and
% it gives unique coordinates to each point in space.
% In the present context, the global system is referenced by uppercase
% letters.
%%%
% *Local system*
%
% The local system xyz is a relative reference. Its origin is located at
% the beginning of each element (initial node) and the local x-axis is
% always in the same direction of the element longitudinal axis.
% Local y-axis and z-axis are in the principal moment of inertia directions
% of the element cross-section.
% The local system is referenced by lowercase letters.
%
%% Analysis Model Types
%
% The LESM program considers only static linear-elastic structural analysis
% of 2D (plane) linear elements models and 3D (spatial) linear elements
% models, which can be of the following types:
%%%
% *Truss analysis*
%
% A truss model is a common form of analysis model, with the following
% assumptions:
%%%
% * Truss elements are bars connected at their ends only, and they are
%   connected by friction-less pins. Therefore, a truss element does not
%   present any secondary bending moment or torsion moment induced by
%   rotation continuity at joints.
% * A truss model is loaded only at nodes.
%   Any load action along an element, such as self weight, is statically
%   transferred as concentrated forces to the element end nodes.
% * Local bending of elements due to internal loads is neglected, when
%   compared to the effect of global load acting on the truss.
% * Therefore, there is only one type of internal force in a truss
%   element: axial internal force, which may be tension or compression.
% * A 2D truss model is considered to be laid in the XY-plane, with only
%   in-plane behavior, that is, there is no displacement transversal to
%   the truss plane.
% * Each node of a 2D truss model has two d.o.f.'s (degrees of
%   freedom): a horizontal displacement in X direction and a vertical
%   displacement in Y direction.
% * Each node of a 3D truss model has three d.o.f.'s (degrees of
%   freedom): displacements in X, Y and Z directions.
%%%
% *Frame analysis*
%
% These are the basic assumptions of a frame model:
%%%
% * Frame elements are usually rigidly connected at joints.
%   However, a frame element might have a hinge (rotation liberation)
%   at an end or hinges at both ends.
% * It is assumed that a hinge in a 2D frame element releases continuity of
%   rotation in the out-of-plane direction, while a hinge in a 3D frame
%   element releases continuity of rotation in all directions.
% * A 2D frame model is considered to be laid in the XY-plane, with only
%   in-plane behavior, that is, there is no displacement transversal to
%   the frame plane.
% * Internal forces at any cross-section of a 2D frame element are:
%   axial force (local x direction), shear force (local y direction),
%   and bending moment (about local z direction).
% * Internal forces at any cross-section of a 3D frame element are:
%   axial force, shear forces (local y direction and local z direction),
%   bending moments (about local y direction and local z direction),
%   and torsion moment (about x direction).
% * Each node of a 2D frame model has three d.o.f.'s: a horizontal
%   displacement in X direction, a vertical displacement in Y
%   direction, and a rotation about the Z-axis.
% * Each node of a 3D frame model has six d.o.f.'s: displacements in
%   X, Y and Z directions, and rotations about X, Y and Z directions.
%%%
% *Grillage analysis*
%
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
%% Structural Element Types
%
% For frame or grillage models, whose elements have bending effects,
% two types of beam elements are considered:
%
% * Navier (Euler-Bernoulli) beam element.
% * Timoshenko beam element
%
% In Euler-Bernoulli flexural behavior, it is assumed that there is no
% shear deformation. As a consequence, bending of a linear structure
% element is such that its cross-section remains plane and normal to the
% element longitudinal axis.
%
% In Timoshenko flexural behavior, shear deformation is considered in an
% approximated manner. Bending of a linear structure element is such that
% its cross-section remains plane but it is not normal to the element
% longitudinal axis.
%
% In truss models, the two types of elements may be used indistinguishably,
% since there is no bending behavior of a truss element, and
% Euler-Bernoulli elements and Timoshenko elements are equivalent for the
% axial behavior.
%
%% Local Axes of an Element
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
%% Load Types
%
% There are four types of loads considered in the LESM program:
%
% * Concentrated nodal force and moment in global axes directions.
% * Uniformely distributed force on elements, spanning its entire
%   length, in local or in global axes directions.
% * Linearly distributed force on elements, spanning its entire length,
%   in local or in global axes directions.
% * Temperature variation along elements local axes (temperature gradient).
%
% In addition, nodal prescribed displacements and rotations may be specified.
%
% It is assumed that there is a single load case.
%
%% Components of Concentrated Nodal Loads
%
% * In 2D truss models, concentrated nodal loads are force components
%   in global X and global Y directions.
% * In 3D truss models, concentrated nodal loads are force components
%   in global X, global Y and global Z directions.
% * In 2D frame models, concentrated nodal loads are force components
%   in global X and global Y directions, and a moment component
%   about global Z direction.
% * In 3D frame models, concentrated nodal loads are force and moment
%   components in global X, global Y and global Z directions.
% * In grillage models, concentrated nodal loads are a force component in
%   global Z direction, and moment components about global X and global Y
%   directions.
%
%% Components of Distributed Loads on Elements
%
% * In 2D truss or 2D frame models, uniformely or linearly distributed
%   forces have two components, which are in global X and Y directions,
%   or in local x and y directions.
% * In 3D truss or 3D frame models, uniformely or linearly distributed
%   forces have three components, which are in global X, Y and Z directions,
%   or in local x, y and z directions.
% * In grillage models, uniformely or linearly distributed forces have only
%   one component, which is in global Z direction.
%
%% Components of Thermal Loads on Elements
%
% * In 2D truss and 2D frame models, thermal loads are specified by a
%   temperature gradient relative to element local y-axis, and the
%   temperature variation on its center of gravity.
% * In 3D truss and 3D frame models, thermal loads are specified by the
%   temperature gradients relative to element local y and z axes, and the
%   temperature variation on its center of gravity.
% * In grillage models, thermal loads are specified by a temperature
%   gradient relative to element local z-axis, and the temperature
%   variation on its center of gravity.
%
% The temperature gradient relative to an element local axis is the
% difference between the temperature variation on its bottom face (face on
% the negative side of local axis) and its upper face (face on the positive
% side of local axis).
%
%% Materials
%
% All materials in LESM are considered to have linear elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in
% all directions.
%
%% Cross-Sections
%
% All cross-sections in LESM are considered to be of a generic type,
% which means that particular shapes are not specified, only their
% geometric properties are provided, such as area, moment of inertia
% and height.
%
%% Object Oriented Classes
%
% This program adopts an Object Oriented Programming (OOP) paradigm, in
% which the following OOP classes are implemented:
%%%
% * <drv.html Drv: analysis driver class>.
% * <anm.html Anm: analysis model class>.
% * <node.html Node: node class>.
% * <elem.html Elem: element class>.
% * <lelem.html Lelem: load element class>.
% * <material.html Material: material class>.
% * <section.html Section: cross-section class>.
% * <print.html Print: print class>.
%
%% Auxiliary Functions and Files
%%%
% * <include_constants.html include_constants: File with global constants>.
% * <readFile.html readFile: Reads structural model data from a neutral-format file>.
% * <saveFile.html SaveFile: Writes a neutral-format file with structural model data>.
%
%% Non-graphical version
clc;
clear;
close all;

% Specify file name (edit only this line to change file name for target model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'Models/Frame2D/Frame2D_1.lsm';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize analysis driver object
drv = Drv();

% Open input file with model information
fid = fopen(fileName,'rt');

% Check for valid input file
if fid > 0
    % Read model information
    fprintf(1,'Pre-processing...\n');
    [vs,print] = readFile(fid,drv);
    
    % Check input file version compatibility
    if vs == 1
        % Process provided data according to the direct stiffness method
        drv.fictRotConstraint(1); % Create ficticious rotation constraints
        status = drv.process();   % Process data
        drv.fictRotConstraint(0); % Remove ficticious rotation constraints

        % Print analysis results
        if status == 1
            print.results();
        end
    else
        fprintf(1,'This file version is not compatible with the program version!\n');
    end
else 
    fprintf(1,'Error opening input file!\n');
end
