%% Drv (Driver) Class
%
% This is a class in the Object Oriented Programming (OOP) paradigm
% that defines an analysis driver object in the <main.html LESM> (Linear
% Elements Structure Model) program.
%
% The analysis driver object can be interpreted as an object that
% represents the structural model as a whole.
% The properties of a driver object are the main global variables of a
% typical structural analysis program.
% The methods of a driver object are the general steps of the direct
% stiffness method that do not depend on the analysis model or element type.
% Therefore, an analysis driver object is responsible for running the
% functions that drives the structural analysis.
% Only one driver object is created during the analysis process.

%% Class definition
% Definition of super-class *Drv* as a <https://www.mathworks.com/help/matlab/ref/handle-class.html *handle*> class.
classdef Drv < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        anm = [];         % handle to an object of the Anm class
        materials = [];   % vector of handles to objects of the Material class
        sections = [];    % vector of handles to objects of the Section class
        nodes = [];       % vector of handles to objects of the Node class
        elems = [];       % vector of handles to objects of the Elem class
        K = [];           % global stiffness matrix
        F = [];           % global forcing vector
        D = [];           % global displacement vector
        ID = [];          % global d.o.f. numbering matrix
        nmat = 0;         % number of materials
        nsec = 0;         % number of cross-sections
        nnp = 0;          % number of nodes
        nel = 0;          % number of elements
        neq = 0;          % number of equations
        neqfree = 0;      % number of equations of free d.o.f.
        neqfixed = 0;     % number of equations of fixed d.o.f.
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function drv = Drv(anm,mat,sec,nodes,elems,K,F,D,ID,nmat,nsec,nnp,nel,neq,neqfr,neqfx)
            if (nargin > 0)
                drv.anm = anm;
                drv.materials = mat;
                drv.sections = sec;
                drv.nodes = nodes;
                drv.elems = elems;
                drv.K = K;
                drv.F = F;
                drv.D = D;
                drv.ID = ID;
                drv.nmat = nmat;
                drv.nsec = nsec;
                drv.nnp = nnp;
                drv.nel = nel;
                drv.neq = neq;
                drv.neqfree = neqfr;
                drv.neqfixed = neqfx;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Creates or removes fictitious rotation constraints on nodes where
        % all incident elements are hinged, to force global stability
        % during the analysis process by avoiding a singular stifness matrix.
        % Input arguments:
        %  fict: flag to check if constraints must be created or removed
        function fictRotConstraint(drv,fict)
            include_constants;
            
            for n = 1:drv.nnp
                % Create fictitious rotation constraints (before analysis process)
                if fict == 1
                    % Calculate total number of incident elements and hinged
                    % elements on current node
                    [tot,hng] = drv.nodes(n).elemsIncidence(drv);
                    
                    % Check if all incident elements are hinged
                    if tot == hng
                         % Insert ficticious constraint
                         % on free rotation d.o.f.'s.
                        if drv.nodes(n).ebc(4) == FREE_DOF
                            drv.nodes(n).ebc(4) = FICTFIXED_DOF;
                        end
                        if drv.nodes(n).ebc(5) == FREE_DOF
                            drv.nodes(n).ebc(5) = FICTFIXED_DOF;
                        end
                        if drv.nodes(n).ebc(6) == FREE_DOF
                            drv.nodes(n).ebc(6) = FICTFIXED_DOF;
                        end
                    end
                    
                % Remove fictitious rotation constraints (after analysis process)
                elseif fict == 0
                    if drv.nodes(n).ebc(4) == FICTFIXED_DOF
                        drv.nodes(n).ebc(4) = FREE_DOF;
                    end
                    if drv.nodes(n).ebc(5) == FICTFIXED_DOF
                        drv.nodes(n).ebc(5) = FREE_DOF;
                    end
                    if drv.nodes(n).ebc(6) == FICTFIXED_DOF
                        drv.nodes(n).ebc(6) = FREE_DOF;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Dimensions and initializes global stiffness matrix, global forcing
        % vector and global displacement vector.
        function dimKFD(drv)
            drv.K = zeros(drv.neq,drv.neq);
            drv.F = zeros(drv.neq,1);
            drv.D = zeros(drv.neq,1);
        end
        
        %------------------------------------------------------------------
        % Assembles global d.o.f (degree of freedom) numbering ID matrix:
        %  ID(k,n) = equation number of d.o.f. k of node n
        %  Free d.o.f.'s have the initial numbering.
        %  countF --> counts free d.o.f.'s. (numbered first)
        %  countS --> counts fixed d.o.f.'s. (numbered later)
        function assembleDOFNum(drv)
            % Initialize equation numbers for free and fixed d.o.f.'s.
            countF = 0;
            countS = drv.neqfree;
            
            % Check if each d.o.f. is free or fixed to increment eqn.
            % number and store it in the ID matrix.
            for n = 1:drv.nnp
                for k = 1:drv.anm.ndof
                    if drv.ID(k,n) == 0
                        countF = countF + 1;
                        drv.ID(k,n) = countF;
                    else
                        countS = countS + 1;
                        drv.ID(k,n) = countS;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element gather vector (gle) that stores element d.o.f.'s
        % equation numbers.
        function assembleGle(drv)
            for e = 1:drv.nel
                % Initialize element gather vector
                drv.elems(e).gle = zeros(drv.elems(e).nen * drv.anm.ndof, 1);
                i = 0;
                
                % Assemble d.o.f.'s of initial node to element gather vector
                for j = 1:drv.anm.ndof
                    i = i + 1;
                    drv.elems(e).gle(i) = drv.ID(j,drv.elems(e).nodes(1).id);
                end
                
                % Assemble d.o.f.'s of final node to element gather vector
                for j = 1:drv.anm.ndof
                    i = i + 1;
                    drv.elems(e).gle(i) = drv.ID(j,drv.elems(e).nodes(2).id);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assembles global stiffness matrix by adding the contribution of
        % the stiffness matrices of all elements.
        function gblMtx(drv)
            for e = 1:drv.nel
                % Compute element stiffness matrix in global system
                keg = drv.elems(e).gblStiffMtx();
                
                % Assemble element stiffness matrix to global matrix
                drv.assembleElemMtx(keg,e);
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element stiffness matrix coefficients (in global system)
        % to the correct positions of the global stiffness matrix.
        % Input arguments:
        %  keg: element stiffness matrix in global system
        %  e: element identification number
        function assembleElemMtx(drv,keg,e)
            % Extract element gather vector
            gle = drv.elems(e).gle;
            
            % Add contribution of element stifness matrix to
            % global stiffness matrix.
            drv.K(gle,gle) = drv.K(gle,gle) + keg;
            
            % This command is a vector indexing, equivalent to:
            %   ndof = drv.anm.ndof;
            %   nen  = drv.elems(e).nen;
            %
            %   for il = 1:nen*ndof
            %       ig = gle(il);
            %       for jl = 1:nen*ndof
            %           jg = gle(jl);
            %           drv.K(ig,jg) = drv.K(ig,jg) + keg(il,jl);
            %       end
            %   end
            %
            % Where:
            %   il is the local (element) row equation number
            %   jl is the local (element) column equation number
            %   ig is the global row equation number
            %   jg is the global column equation number
        end
        
        %------------------------------------------------------------------
        % Adds element equivalent nodal loads (ENL), from distributed loads
        % and thermal loads, to global forcing vector.
        function elemLoads(drv)
            for e = 1:drv.nel
                % Add distributed load as ENL to global forcing vector
                feg = drv.elems(e).load.gblDistribLoadENL();
                drv.assembleENL(feg,e);
                
                % Add thermal load as ENL to global forcing vector
                feg = drv.elems(e).load.gblThermalLoadENL();
                drv.assembleENL(feg,e);
            end
        end
        
        %------------------------------------------------------------------
        % Assembles element equivalent nodal load vector (in global system)
        % to any term of the global forcing vector, including the terms that
        % correspond to constrained d.o.f.'s.
        % Input arguments:
        %  feg: element equivalent nodal load vector in global system
        %  e: element identification number
        function assembleENL(drv,feg,e)
            % Extract element gather vector
            gle = drv.elems(e).gle;
            
            % Add contribution of element equivalent nodal load vector to
            % global forcing vector
            drv.F(gle) = drv.F(gle) + feg;
            
            % This command is a vector indexing, equivalent to:
            %   ndof = drv.anm.ndof;
            %   nen  = drv.elems(e)nen;
            %
            %   for il = 1:nen*ndof
            %       ig = gle(il);
            %       drv.F(ig) = drv.F(ig) + feg(il);
            %   end
            %
            % Where:
            %   il is the local (element) row equation number
            %   ig is the global row equation number
        end
        
        %------------------------------------------------------------------
        % Partitions and solves the system of equilibrium equations.
        % Partitions the coefficient matrix K, forcing vector F and unknown 
        % vector D:
        %  f -> free d.o.f.'s (numbered first) - natural B.C. (unknown)
        %  s -> constrainted by support d.o.f.'s (numbered later) - essential B.C. (known)
        % 
        % Partitioned equilibrium system:
        % [ Kff Kfs ] * [ Df ] = [ Ff ]
        % [ Ksf Kss ]   [ Ds ] = [ Fs ]
        %
        % Output:
        %  status: flag for stable free-free global matrix (0 = unstable, 1 = stable)
        function status = solveEqnSystem(drv)
            % Assemble free-free global matrix
            Kff = drv.K(1:drv.neqfree,1:drv.neqfree);
            
            % Check for stable free-free global matrix by verifying its
            % reciprocal condition number (a very low reciprocal condition
            % number indicates that the matrix is badly conditioned and
            % may be singular)
            status = 1;
            if (rcond(Kff) < 10e-12)
                status = 0;
                return;
            end
            
            % Partition system of equations
            Kfs = drv.K(1:drv.neqfree,drv.neqfree+1:drv.neq);
            Ksf = drv.K(drv.neqfree+1:drv.neq,1:drv.neqfree);
            Kss = drv.K(drv.neqfree+1:drv.neq,drv.neqfree+1:drv.neq);
            Ff  = drv.F(1:drv.neqfree);
            Ds  = drv.D(drv.neqfree+1:drv.neq);
            Fs  = drv.F(drv.neqfree+1:drv.neq);
            
            % Solve for Df
            Df = Kff \ (Ff - Kfs * Ds);
            
            % Reconstruct the global unknown vector D
            drv.D = [ Df
                      Ds ];
            
            % Recover forcing unknown values (reactions) at essential B.C.
            % It is assumed that the Fs vector currently stores combined 
            % nodal loads applied directly to fixed d.o.f's.
            % Superimpose computed reaction values to combined nodal loads,
            % with inversed direction, that were applied directly to fixed 
            % d.o.f.'s.
            Fs = -Fs + Ksf * Df + Kss * Ds;
            
            % Reconstruct the global forcing vector
            drv.F = [ Ff
                      Fs ];
        end
        
        %------------------------------------------------------------------
        % Computes element internal forces on element ends.
        function elemIntForce(drv)
            for e = 1:drv.nel
                % Initialize element internal force arrays with null values
                drv.anm.initIntForce(drv.elems(e));
                
                % Compute element internal forces from global analysis
                % (from nodal displacements and rotations)
                fel_gblAnl = drv.elems(e).gblAnlIntForce(drv);
                
                % Get element internal forces from local analysis
                % (fixed end forces of distributed loads and thermal loads)
                fel_lclAnl = drv.elems(e).fel_distribLoad + drv.elems(e).fel_thermalLoad;
                
                % Add contribution of element internal forces from global
                % analysis and from local analysis
                fel = fel_gblAnl + fel_lclAnl;
                
                % Assemble element internal force arrays
                drv.anm.assembleIntForce(drv.elems(e),fel);
            end
        end
        
        %------------------------------------------------------------------
        % Computes element internal displacements in several cross-section
        % positions.
        function elemIntDispl(drv)
            for e = 1:drv.nel
                % Initialize element internal displacements array with
                % null values
                drv.anm.initIntDispl(drv.elems(e));
                
                % Compute displacements in several cross-section positions
                L = drv.elems(e).length;
                step = L / (size(drv.elems(e).intDispl,2) - 1);
                j = 1;
                for x = 0:step:L
                    % Compute element internal displacements from global
                    % analysis (from nodal displacements and rotations)
                    del_gblAnl = drv.elems(e).gblAnlIntDispl(drv,x);
                    
                    % Compute element internal displacements from local
                    % analysis (internal displacements from distributed
                    % loads and thermal loads)
                    del_lclAnl = drv.anm.lclAnlIntDispl(drv.elems(e),x);
                    
                    % Add contribution of element internal displacements
                    % from global analysis and from local analysis
                    del = del_gblAnl + del_lclAnl;
                    
                    % Assemble element internal displacements array
                    drv.elems(e).intDispl(:,j) = del;
                    j = j + 1;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Processes current model data according to the direct stiffness
        % method, i.e., assembles global equilibrium system of equations,
        % calculates nodal displacements and rotations, and computes
        % support reactions and elements internal forces.
        % This method also prints the processing stages and results in 
        % the default output (MATLAB command window).
        % Output:
        %  status: flag for stable free-free global matrix (0 = unstable, 1 = stable)
        function status = process(drv)
            fprintf(1,'Preparing analysis data...\n');
            
            % Dimension and initialize global matrix and vectors
            drv.dimKFD();
            
            % Generate global d.o.f. numbering matrix
            drv.anm.setupDOFNum(drv);
            
            % Assemble global d.o.f. numbering matrix
            drv.assembleDOFNum();
            
            % Assemble element gather vectors
            drv.assembleGle();
            
            % Store prescribed displacements in global displacement vector
            drv.anm.setupPrescDispl(drv);
            
            % Assemble global stiffness matrix
            fprintf(1,'Assembling stiffness matrix...\n');
            drv.gblMtx();
            
            % Assemble global forcing vector
            fprintf(1,'Assembling forcing vector...\n');
            drv.anm.nodalLoads(drv); % Initialize forcing vector with nodal loads
            drv.elemLoads();         % Add element equivalent nodal loads to forcing vector
            
            % Partition and solve the system of equations
            fprintf(1,'Solving system of equations...\n');
            status = drv.solveEqnSystem();
            
            % Check structural model stability
            if (status == 1)
                % Compute elements internal forces
                fprintf(1,'Computing element internal forces...\n');
                drv.elemIntForce();
                
                % Compute elements internal displacements
                fprintf(1,'Computing element internal displacements...\n');
                drv.elemIntDispl();
            else
                fprintf(1,'Unstable structure or invalid input data...\n');
            end
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of the entire model.
        function clean(drv)
            % Clean data structure of the Anm object
            if drv.anm ~= 0
                drv.anm.clean();
            end
            
            % Clean data structure of all Material objects
            for m = 1:drv.nmat
                drv.materials(m).clean();
            end
            
            % Clean data structure of all Section objects
            for s = 1:drv.nsec
                drv.sections(s).clean();
            end
            
            % Clean data structure of all Node objects
            for n = 1:drv.nnp
                drv.nodes(n).clean();
            end
            
            % Clean data structure of all Element objects and Lelem objects
            for e = 1:drv.nel
                drv.elems(e).load.clean();
                drv.elems(e).clean();
            end
            
            % Clean data structure of the Drv object
            drv.anm = [];
            drv.materials = [];
            drv.sections = [];
            drv.nodes = [];
            drv.elems = [];
            drv.K = [];
            drv.F = [];
            drv.D = [];
            drv.ID = [];
            drv.nmat = 0;
            drv.nsec = 0;
            drv.nnp = 0;
            drv.nel = 0;
            drv.neq = 0;
            drv.neqfree = 0;
            drv.neqfixed = 0;
        end
    end
end