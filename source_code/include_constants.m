%% Global Constants
% This file contains variables that have a global scope and store flags to
% help understanding the <main.html LESM> program.
% It is necessary to include this file in every function or method that
% uses a global constant.
%
%% Types of analysis models
global TRUSS2D_ANALYSIS        % 2D truss analysis
global FRAME2D_ANALYSIS        % 2D frame analysis
global GRILLAGE_ANALYSIS       % grillage analysis
global TRUSS3D_ANALYSIS        % 3D truss analysis
global FRAME3D_ANALYSIS        % 3D frame analysis

TRUSS2D_ANALYSIS    = 0;
FRAME2D_ANALYSIS    = 1;
GRILLAGE_ANALYSIS   = 2;
TRUSS3D_ANALYSIS    = 3;
FRAME3D_ANALYSIS    = 4;

%% Types of elements
global MEMBER_NAVIER           % Navier (Euler-Bernoulli) beam element
global MEMBER_TIMOSHENKO       % Timoshenko beam element

MEMBER_NAVIER       = 0;
MEMBER_TIMOSHENKO   = 1;

%% Types of continuity conditions
global HINGED_END              % hinged element end
global CONTINUOUS_END          % continuous element end

HINGED_END          = 0;
CONTINUOUS_END      = 1;

%% Types of load directions
global GLOBAL_LOAD             % element load applied in global system
global LOCAL_LOAD              % element load applied in local system

GLOBAL_LOAD         = 0;
LOCAL_LOAD          = 1;

%% Types of degree of freedom constraints
global FREE_DOF                % free degree of freedom
global FIXED_DOF               % fixed degree of freedom
global FICTFIXED_DOF           % fictitious rotation constraints

FREE_DOF            =  0;
FIXED_DOF           =  1;
FICTFIXED_DOF       =  2;
