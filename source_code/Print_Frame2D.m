%% Print_Frame2D Class
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <print.html *Print*> in the <main.html LESM> (Linear
% Elements Structure Model) program. This sub-class implements abstract
% methods, declared in super-class *Print*, that deal with 2D frame
% analysis model of linear elements.
%
%% Class definition
% Definition of sub-class *Print_Frame2D* derived from super-class <print.html Print>.
classdef Print_Frame2D < Print
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function print = Print_Frame2D(drv)
            print = print@Print(1);
            
            if (nargin > 0)
                print.drv = drv;
            end
        end
    end

    %% Public methods
    % Implementation of the abstract methods declared in super-class <print.html *Print*>.
    methods
        %------------------------------------------------------------------
        % Prints analyis results.
        function results(print)
            print.header();
            fprintf(print.txt, '\n\n\n____________ M O D E L  I N F O R M A T I O N ____________\n');
            print.analysisLabel();
            [n_nl,n_pd,n_ul,n_ll,n_tv] = print.modelDescrip();
            print.material();
            print.section();
            print.nodalCoords();
            print.nodalSupport();
            print.nodalLoads(n_nl);
            print.nodalPrescDisp(n_pd);
            print.elements();
            print.unifElementLoads(n_ul);
            print.linearElementLoads(n_ll);
            print.temperatureVariation(n_tv);
            fprintf(print.txt, '\n\n\n\n_____________ A N A L Y S I S  R E S U L T S _____________\n');
            print.nodalDisplRot();
            print.reactions();
            print.intForces();
            print.elemDispl();
        end

        %------------------------------------------------------------------
        % Prints analyis model type.
        function analysisLabel(print)
            fprintf(print.txt, '\n\n----------------------------\n');
            fprintf(print.txt, 'ANALYSIS MODEL: ');
            fprintf(print.txt, '2D FRAME\n');
            fprintf(print.txt, '----------------------------\n');
        end
        
        %------------------------------------------------------------------
        % Prints model description.
        % Output:
        %  n_nl: number of nodes with applied loads
        %  n_pd: number of nodes with prescribed displacement
        %  n_ul: number of elements with uniformly distributed loads
        %  n_ll: number of elements with linearly distributed loads
        %  n_tv: number of elements with temperature variation
        function [n_nl,n_pd,n_ul,n_ll,n_tv] = modelDescrip(print)
            include_constants;
            
            % Initialize variables
            neqfixed = 0;
            n_nl = 0;
            n_pd = 0;
            n_ul  = 0;
            n_ll = 0;
            n_tv = 0;
            
            % Loop over all nodes
            for n = 1:print.drv.nnp
                % Increment number of fixed d.o.f.
                if print.drv.nodes(n).ebc(1) == FIXED_DOF
                    neqfixed = neqfixed + 1;
                end
                if print.drv.nodes(n).ebc(2) == FIXED_DOF
                    neqfixed = neqfixed + 1;
                end
                if print.drv.nodes(n).ebc(6) == FIXED_DOF
                    neqfixed = neqfixed + 1;
                end
                
                % Increment number of nodes with applied loads
                if ~isempty(print.drv.nodes(n).nodalLoad)
                    n_nl = n_nl + 1;
                end
                
                % Increment number of nodes with prescribed displacement
                if ~isempty(print.drv.nodes(n).prescDispl)
                    n_pd = n_pd + 1;
                end
            end
            neqfree = print.drv.neq - neqfixed; % number of free d.o.f.
            
            % Loop over all elements
            for e = 1:print.drv.nel
                % Increment number of elements with uniformly distributed loads
                if ~isempty(print.drv.elems(e).load.uniformGbl)
                    n_ul = n_ul + 1;
                end
                
                % Increment number of elements with linearly distributed loads
                if ~isempty(print.drv.elems(e).load.linearGbl)
                    n_ll = n_ll + 1;
                end
                
                % Increment number of elements with temperature variation
                if (print.drv.elems(e).load.tempVar_X ~= 0) || ...
                   (print.drv.elems(e).load.tempVar_Y ~= 0) || ...
                   (print.drv.elems(e).load.tempVar_Z ~= 0)
                    n_tv = n_tv + 1;
                end
            end

            fprintf(print.txt, '\n\n----------------------------------------\n' );
            fprintf(print.txt, 'M O D E L  D E S C R I P T I O N:\n' );
            fprintf(print.txt, '----------------------------------------\n');
            fprintf(print.txt, 'NUMBER OF NODES.....................:%4d\n', print.drv.nnp);
            fprintf(print.txt, 'NUMBER OF ELEMENTS .................:%4d\n', print.drv.nel);
            fprintf(print.txt, 'NUMBER OF DEGREES OF FREEDOM........:%4d\n', print.drv.neq);
            fprintf(print.txt, 'NUMBER OF FREE DEGREES OF FREEDOM...:%4d\n', neqfree);
            fprintf(print.txt, 'NUMBER OF FIXED DEGREES OF FREEDOM..:%4d\n', neqfixed);
            fprintf(print.txt, 'NUMBER OF MATERIALS.................:%4d\n', print.drv.nmat);
            fprintf(print.txt, 'NUMBER OF CROSS-SECTIONS............:%4d\n', print.drv.nsec);
            fprintf(print.txt, 'NUMBER OF NODES W/ APPLIED LOADS....:%4d\n', n_nl);
            fprintf(print.txt, 'NUMBER OF NODES W/ PRESC. DISPL.....:%4d\n', n_pd);
            fprintf(print.txt, 'NUMBER OF ELEMENTS W/ UNIFORM LOAD..:%4d\n', n_ul);
            fprintf(print.txt, 'NUMBER OF ELEMENTS W/ LINEAR LOAD...:%4d\n', n_ll);
            fprintf(print.txt, 'NUMBER OF ELEMENTS W/ TEMP. VAR.....:%4d\n', n_tv);
            fprintf(print.txt, '----------------------------------------\n');
        end
        
        %------------------------------------------------------------------
        % Prints cross-section properties.
        function section(print)
            fprintf(print.txt, '\n\n----------------------------------------------\n');
            fprintf(print.txt, 'C R O S S  S E C T I O N  P R O P E R T I E S\n');
            fprintf(print.txt, '----------------------------------------------\n');

            if print.drv.nsec > 0
                fprintf(print.txt, ' SECTION     FULL AREA [cm²]     SHEAR AREA [cm²]     INERTIA [cm4]      HEIGHT [cm]\n');

                for s = 1:print.drv.nsec
                    fprintf(print.txt, '%4d   %15.2f   %17.2f   %18.2f   %15.2f\n', s,...
                            1e4*print.drv.sections(s).area_x,...
                            1e4*print.drv.sections(s).area_y,...
                            1e8*print.drv.sections(s).inertia_z,...
                            1e2*print.drv.sections(s).height_y);
                end
            else
                fprintf(print.txt, ' NO CROSS-SECTION\n');
            end
        end

        %------------------------------------------------------------------
        % Prints nodal support conditions.
        function nodalSupport(print)
            include_constants;
            
            fprintf(print.txt, '\n\n------------------------------\n');
            fprintf(print.txt, 'N O D A L  R E S T R A I N T S \n');
            fprintf(print.txt, '------------------------------\n');
            fprintf(print.txt, ' NODE   DISPL X   DISPL Y   ROTAT Z\n');

            for n = 1:print.drv.nnp
                if print.drv.nodes(n).ebc(1) == FIXED_DOF
                    node_restr1 = 'FIXED';
                else
                    node_restr1 = 'FREE';
                end

                if print.drv.nodes(n).ebc(2) == FIXED_DOF
                    node_restr2 = 'FIXED';
                else
                    node_restr2 = 'FREE';
                end

                if print.drv.nodes(n).ebc(6) == FIXED_DOF
                    node_restr3 = 'FIXED';
                else
                    node_restr3 = 'FREE';
                end

                fprintf(print.txt, '%4d     %5s     %5s     %5s\n', n,...
                        node_restr1,...
                        node_restr2,...
                        node_restr3);
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal loads.
        % Input arguments:
        %  n_nl: number of nodes with applied loads
        function nodalLoads(print,n_nl)
            fprintf(print.txt, '\n\n--------------------\n');
            fprintf(print.txt, 'N O D A L  L O A D S \n');
            fprintf(print.txt, '--------------------\n');
            
            if n_nl ~= 0
                fprintf(print.txt, ' NODE     FX [kN]         FY [kN]         MZ [kNm]\n');

                for n = 1:print.drv.nnp
                    if ~isempty(print.drv.nodes(n).nodalLoad)
                        fprintf(print.txt, '%4d    %9.3f    %12.3f    %12.3f\n', n,...
                                print.drv.nodes(n).nodalLoad(1),...
                                print.drv.nodes(n).nodalLoad(2),...
                                print.drv.nodes(n).nodalLoad(6));
                    end
                end
            else
                fprintf(print.txt, ' NO NODAL LOAD\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints nodal prescribed displacements.
        % Input arguments:
        %  n_pd: number of nodes with prescribed displacement
        function nodalPrescDisp(print,n_pd)
            fprintf(print.txt, '\n\n--------------------------------\n');
            fprintf(print.txt, 'N O D A L  P R E S C.  D I S P L. \n');
            fprintf(print.txt, '--------------------------------\n');
            
            if n_pd ~= 0
                fprintf(print.txt, ' NODE     DX [mm]          DY [mm]          RZ [rad]\n');

                for n = 1:print.drv.nnp
                    if ~isempty(print.drv.nodes(n).prescDispl)
                        fprintf(print.txt, '%4d   %8.1f   %14.1f   %16.3f\n', n,...
                                1e3*print.drv.nodes(n).prescDispl(1),...
                                1e3*print.drv.nodes(n).prescDispl(2),...
                                print.drv.nodes(n).prescDispl(6) );
                    end
                end
            else
                fprintf(print.txt, ' NO PRESCRIBED DISPLACEMENT\n');
            end
        end
        
        %------------------------------------------------------------------
        % Prints elements information.
        function elements(print)
            include_constants;

            fprintf(print.txt, '\n\n---------------\n');
            fprintf(print.txt, 'E L E M E N T S\n');
            fprintf(print.txt, '---------------\n');

            if print.drv.nel > 0
            
                for e = 1:print.drv.nel
                    element_type = print.drv.elems(e).type;

                    switch element_type
                        case MEMBER_NAVIER
                            type = 'Navier element';
                        case MEMBER_TIMOSHENKO
                            type = 'Timoshenko element';
                    end

                    if print.drv.elems(e).hingei == HINGED_END
                        hingei = 'yes';
                    else
                        hingei = 'no ';
                    end

                    if print.drv.elems(e).hingef == HINGED_END
                        hingef = 'yes';
                    else
                        hingef = 'no ';
                    end

                    mat = print.drv.elems(e).material.id;
                    sec = print.drv.elems(e).section.id;
                    nodei = print.drv.elems(e).nodes(1).id;
                    nodef = print.drv.elems(e).nodes(2).id;
                    L = print.drv.elems(e).length;
                    v = [print.drv.elems(e).vz(1), print.drv.elems(e).vz(2), print.drv.elems(e).vz(3)];
                    v = v / norm(v);
                    
                    if strcmp(type,'Navier element') && e == 1
                    fprintf(print.txt, ' ELEMENT      TYPE       MAT  SEC  HINGEi HINGEf  NODEi NODEf  LENGTH [m]   vz_X     vz_Y     vz_Z\n');
                    elseif strcmp(type,'Timoshenko element') && e == 1
                    fprintf(print.txt, ' ELEMENT      TYPE          MAT  SEC  HINGEi HINGEf  NODEi NODEf  LENGTH [m]   vz_X     vz_Y     vz_Z\n');    
                    end
                    
                    fprintf(print.txt, '%5d   %15s %4d %4d   %3s    %3s    %3d  %4d  %9.3f     %6.3f   %6.3f  %6.3f\n', ...
                            e, type, mat, sec, hingei, hingef, nodei, nodef, L, v(1), v(2), v(3));
                end
            else
                fprintf(print.txt, ' NO ELEMENT\n');
            end
        end

        %------------------------------------------------------------------
        % Prints uniformly distributed loads information.
        % Input arguments:
        %  n_ul: number of elements with uniformly distributed loads
        function unifElementLoads(print,n_ul)
            include_constants;
            
            fprintf(print.txt, '\n\n----------------------------------------\n');
            fprintf(print.txt, 'U N I F O R M  E L E M E N T  L O A D S\n');
            fprintf(print.txt, '----------------------------------------\n');
            
            if n_ul ~= 0
                fprintf(print.txt, ' ELEMENT  DIRECTION     QX [kN/m]     QY [kN/m]\n');
                
                for e = 1:print.drv.nel
                    if ~isempty(print.drv.elems(e).load.uniformGbl)
                        if print.drv.elems(e).load.uniformDir == GLOBAL_LOAD
                            dir = 'GLOBAL';
                            qx = print.drv.elems(e).load.uniformGbl(1);
                            qy = print.drv.elems(e).load.uniformGbl(2);
                        else
                            dir = 'LOCAL';
                            qx = print.drv.elems(e).load.uniformLcl(1);
                            qy = print.drv.elems(e).load.uniformLcl(2);
                        end

                        fprintf(print.txt, '%5d   %9s %13.3f %13.3f\n',e,dir,qx,qy);
                    end
                end
            else
                fprintf(print.txt, ' NO UNIFORM LOAD\n');
            end
        end

        %------------------------------------------------------------------
        % Prints linearly distributed loads information.
        % Input arguments:
        %  n_ll: number of elements with linearly distributed loads
        function linearElementLoads(print,n_ll)
            include_constants;
            
            fprintf(print.txt, '\n\n--------------------------------------\n');
            fprintf(print.txt, 'L I N E A R  E L E M E N T  L O A D S\n');
            fprintf(print.txt, '--------------------------------------\n');
            
            if n_ll ~= 0
                fprintf(print.txt, ' ELEMENT  DIRECTION     QXi [kN/m]     QXf [kN/m]     QYi [kN/m]     QYf [kN/m]\n');

                for e = 1:print.drv.nel
                    if ~isempty(print.drv.elems(e).load.linearGbl)
                        if print.drv.elems(e).load.linearDir == GLOBAL_LOAD
                            dir = 'GLOBAL';
                            qxi = print.drv.elems(e).load.linearGbl(1);
                            qyi = print.drv.elems(e).load.linearGbl(2);
                            qxf = print.drv.elems(e).load.linearGbl(4);
                            qyf = print.drv.elems(e).load.linearGbl(5);
                        else
                            dir = 'LOCAL ';
                            qxi = print.drv.elems(e).load.linearLcl(1);
                            qyi = print.drv.elems(e).load.linearLcl(2);
                            qxf = print.drv.elems(e).load.linearLcl(4);
                            qyf = print.drv.elems(e).load.linearLcl(5);
                        end

                        fprintf(print.txt, '%5d   %9s %14.3f %14.3f %14.3f %14.3f\n',...
                                e, dir, qxi, qxf, qyi, qyf);
                    end
                end
            else
                fprintf(print.txt, ' NO LINEAR LOAD\n');
            end
        end

        %------------------------------------------------------------------
        % Prints thermal loads information.
        % Input arguments:
        %  n_tv: number of elements with temperature variation
        function temperatureVariation(print,n_tv)
            fprintf(print.txt, '\n\n-----------------------------------------\n');
            fprintf(print.txt, 'T E M P E R A T U R E  V A R I A T I O N\n');
            fprintf(print.txt, '-----------------------------------------\n');
            
            if n_tv ~= 0
                fprintf(print.txt, ' ELEMENT     dTx [°C]     dTy [°C]\n');

                for e = 1:print.drv.nel
                    if (print.drv.elems(e).load.tempVar_X ~= 0) || ...
                       (print.drv.elems(e).load.tempVar_Y ~= 0)
                        dtx = print.drv.elems(e).load.tempVar_X;
                        dty = print.drv.elems(e).load.tempVar_Y;

                        fprintf(print.txt, '%5d   %12.3f %12.3f\n', e, dtx, dty);
                    end
                end
            else
                fprintf(print.txt, ' NO TEMPERATURE VARIATION\n');
            end
        end

        %------------------------------------------------------------------
        % Prints results of nodal displacement/rotation.
        function nodalDisplRot(print)
            fprintf(print.txt, '\n\n---------------------------------------------------------\n');
            fprintf(print.txt, 'N O D A L  D I S P L A C E M E N T S / R O T A T I O N S \n');
            fprintf(print.txt, '---------------------------------------------------------\n');
            fprintf(print.txt, ' NODE       DISPL X [mm]       DISPL Y [mm]       ROTAT Z [rad]\n');

            for n = 1:print.drv.nnp
                dx = 1e3*print.drv.D(print.drv.ID(1,n));
                dy = 1e3*print.drv.D(print.drv.ID(2,n));
                rz = print.drv.D(print.drv.ID(3,n));
                fprintf(print.txt, '%4d     %11.2f     %14.2f     %17.5f\n', n, dx, dy, rz);
            end
        end

        %------------------------------------------------------------------
        % Prints results of support reactions.
        function reactions(print)
            fprintf(print.txt, '\n\n---------------------------------\n');
            fprintf(print.txt, 'S U P P O R T  R E A C T I O N S\n');
            fprintf(print.txt, '---------------------------------\n');
            fprintf(print.txt, ' NODE       FORCE X [kN]     FORCE Y [kN]     MOMENT Z [kNm]\n');

            for n = 1:print.drv.nnp
                if( (print.drv.ID(1,n) > print.drv.neqfree) || ...
                    (print.drv.ID(2,n) > print.drv.neqfree) || ...
                    (print.drv.ID(3,n) > print.drv.neqfree) )

                    if print.drv.ID(1,n) > print.drv.neqfree 
                        node_reaction1 = print.drv.F(print.drv.ID(1,n));
                    else
                        node_reaction1 = 0.0;
                    end

                    if print.drv.ID(2,n) > print.drv.neqfree 
                        node_reaction2 = print.drv.F(print.drv.ID(2,n));
                    else
                        node_reaction2 = 0.0;
                    end

                    if print.drv.ID(3,n) > print.drv.neqfree
                        node_reaction3 = print.drv.F(print.drv.ID(3,n));
                    else
                        node_reaction3 = 0.0;
                    end

                    fprintf(print.txt, '%4d     %12.3f   %14.3f   %15.3f\n', n, ...
                            node_reaction1, node_reaction2, node_reaction3);
                end
            end
        end

        %------------------------------------------------------------------
        % Prints results of internal forces at element nodes.
        function intForces(print)
            fprintf(print.txt, '\n\n------------------------------------------------------------\n');
            fprintf(print.txt, 'I N T E R N A L  F O R C E S  A T  E L E M E N T  N O D E S\n');
            fprintf(print.txt, '------------------------------------------------------------\n');
            fprintf(print.txt, ' ELEMENT   AXIAL FORCE [kN]       SHEAR FORCE [kN]      BENDING MOMENT [kNm]\n');
            fprintf(print.txt, '           Nodei      Nodef        Nodei     Nodef        Nodei      Nodef\n');

            for e = 1:print.drv.nel
                fprintf(print.txt, '%4d  %10.3f %10.3f   %10.3f %9.3f   %10.3f %10.3f\n', e, ...
                        print.drv.elems(e).axial_force(1),...
                        print.drv.elems(e).axial_force(2), ...
                        print.drv.elems(e).shear_force_Y(1),...
                        print.drv.elems(e).shear_force_Y(2), ...
                        print.drv.elems(e).bending_moment_Z(1),...
                        print.drv.elems(e).bending_moment_Z(2));
            end
        end
        
        %------------------------------------------------------------------
        % Prints results of elements internal displacements.
        function elemDispl(print)
            fprintf(print.txt, '\n\n-----------------------------------------------------------------------------------------\n');
            fprintf(print.txt, 'E L E M E N T S  I N T E R N A L  D I S P L A C E M E N T S  I N  L O C A L  S Y S T E M\n' );
            fprintf(print.txt, '-----------------------------------------------------------------------------------------\n');
            fprintf(print.txt, 'Axial and transversal displacements in 10 cross-sections from x = 0 to x = L\n\n');

            for e = 1:print.drv.nel
                L = print.drv.elems(e).length;
                l = zeros(10,1);
                l(1) = 0;
                for i = 2:10
                    l(i) = l(i-1) + L/9;
                end
                
                du = zeros(10,1);
                dv = zeros(10,1);
                j = 1;
                for x = 0:L/9:L
                    del_gblAnl = 1000*print.drv.elems(e).gblAnlIntDispl(print.drv,x);
                    del_lclAnl = 1000*print.drv.anm.lclAnlIntDispl(print.drv.elems(e),x);
                    d = del_gblAnl + del_lclAnl;
                    du(j) = d(1);
                    dv(j) = d(2);
                    j = j + 1;
                end
                
                fprintf(print.txt, ' ELEM %d\n', e);
                
                fprintf(print.txt, '   X [m]     %6.3f     %6.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f   %8.3f\n',...
                                                 l(1),     l(2),     l(3),    l(4),    l(5),    l(6),    l(7),    l(8),    l(9),   l(10));
                
                fprintf(print.txt, ' du [mm]  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f\n',...
                                               du(1), du(2),  du(3),  du(4),   du(5), du(6),   du(7),  du(8), du(9),   du(10));
                                           
                fprintf(print.txt, ' dv [mm]  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f  %+9.3f\n\n',...
                                               dv(1),  dv(2),  dv(3),  dv(4),  dv(5),  dv(6),  dv(7),  dv(8),  dv(9),  dv(10));
            end
        end
    end
end