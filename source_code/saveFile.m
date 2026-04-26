%% Save File Function
%
% This is an auxiliary function of the <main.html LESM>
% (Linear Elements Structure Model) program that writes all information
% about a structural model in a neutral-format file with the _.lsm_
% extension.
%
%% Function Code
function saveFile(drv,lsm)
include_constants;

fprintf(lsm, 'NEUTRAL-FORMAT FILE\n');
fprintf(lsm, 'This file stores all information about a structure model that can be read and\n');
fprintf(lsm, 'loaded by the LESM (Linear Elements Structure Model) program.\n');
fprintf(lsm, 'To modify the model, edit only the data below the %%TAGS\n\n');

%--------------------------------------------------------------------------
% Version header
fprintf(lsm, '-----------------------------------------------------------------------\n');
fprintf(lsm, 'Specify version number: ''X.XX''\n');
fprintf(lsm, '%%HEADER.VERSION\n');
fprintf(lsm, '1.1');
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Analysis header
fprintf(lsm, '-----------------------------------------------------------------------\n');
fprintf(lsm, 'Specify analysis model type: ''TRUSS2D'', ''FRAME2D'', ''Grillage'', ''TRUSS3D'' OR ''FRAME3D''\n');
fprintf(lsm, '%%HEADER.ANALYSIS\n');
mdata = guidata(findobj('Tag','GUI_Main'));
anm = get(mdata.popupmenu_Anm,'Value') - 1;
if anm == TRUSS2D_ANALYSIS
    fprintf(lsm, '''TRUSS2D''');
elseif anm == FRAME2D_ANALYSIS
    fprintf(lsm, '''FRAME2D''');
elseif anm == GRILLAGE_ANALYSIS
    fprintf(lsm, '''GRILLAGE''');
elseif anm == TRUSS3D_ANALYSIS
    fprintf(lsm, '''TRUSS3D''');
elseif anm == FRAME3D_ANALYSIS
    fprintf(lsm, '''FRAME3D''');
end
fprintf(lsm, '\n\n');

%--------------------------------------------------------------------------
% Nodal coordinates
if drv.nnp > 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Provide nodal coordinates\n');
    fprintf(lsm, 'First line: Total number of nodes\n');
    fprintf(lsm, 'Following lines: Node ID, coord_X [m], coord_Y [m], coord_Z [m]\n');
    fprintf(lsm, '%%NODE.COORD\n');
    fprintf(lsm, '%d\n', drv.nnp);
    for n = 1:drv.nnp
        fprintf(lsm, '%d     %d  %d  %d\n', n, drv.nodes(n).coord(1),...
                                               drv.nodes(n).coord(2),...
                                               drv.nodes(n).coord(3));
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Support conditions
if drv.nnp > 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Specify essential boundary conditions (support conditions)\n');
    fprintf(lsm, 'First line: Total number of nodes\n');
    fprintf(lsm, 'Following lines: Node ID, ebc_dX, ebc_dY, ebc_dZ, ebc_rX, ebc_rY, ebc_rZ\n');
    fprintf(lsm, 'ebc = 0 --> Free degree-of-freedom\n');
    fprintf(lsm, 'ebc = 1 --> Fixed degree-of-freedom\n');
    fprintf(lsm, '%%NODE.SUPPORT\n');
    fprintf(lsm, '%d\n', drv.nnp);
    for n = 1:drv.nnp
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d\n', n,...
                drv.nodes(n).ebc(1), drv.nodes(n).ebc(2),...
                drv.nodes(n).ebc(3), drv.nodes(n).ebc(4),...
                drv.nodes(n).ebc(5), drv.nodes(n).ebc(6));
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Material properties
if drv.nmat > 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Provide materials\n');
    fprintf(lsm, 'First line: Total number of materials\n');
    fprintf(lsm, 'Following lines: Material ID, Elasticity [MPa], Poisson Ratio, Thermal exp. coeff. [/oC]\n');
    fprintf(lsm, '%%MATERIAL.ISOTROPIC\n');
    fprintf(lsm, '%d\n', drv.nmat);
    for m = 1:drv.nmat
        fprintf(lsm, '%d     %d  %d  %d\n', m,...
                1e-3 * drv.materials(m).elasticity,...
                drv.materials(m).poisson,...
                drv.materials(m).thermExp);
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Cross-section properties
if drv.nsec > 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Provide cross-sections\n');
    fprintf(lsm, 'First line: Total number of cross-sections\n');
    fprintf(lsm, 'Following lines: Cross-section ID, Area_x [cm2], Area_y [cm2], Area_z [cm2]\n');
    fprintf(lsm, '                 Inertia_x [cm4], Inertia_y [cm4], Inertia_z [cm4]\n');
    fprintf(lsm, '                 Height_y [cm], Height_z [cm]\n');
    fprintf(lsm, '%%SECTION.PROPERTY\n');
    fprintf(lsm, '%d\n', drv.nsec);
    for s = 1:drv.nsec
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d  %d  %d\n', s,...
                1e4*drv.sections(s).area_x,...
                1e4*drv.sections(s).area_y,...
                1e4*drv.sections(s).area_z,...
                1e8*drv.sections(s).inertia_x,...
                1e8*drv.sections(s).inertia_y,...
                1e8*drv.sections(s).inertia_z,...
                1e2*drv.sections(s).height_y,...
                1e2*drv.sections(s).height_z);
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Beam end liberation
fprintf(lsm, '-----------------------------------------------------------------------\n');
fprintf(lsm, 'Flag specification for different beam end liberation conditions\n');
fprintf(lsm, '1 - Continuous end / Continuous end\n');
fprintf(lsm, '2 - Hinged end / Continuous end\n');
fprintf(lsm, '3 - Continuous end / Hinged end\n');
fprintf(lsm, '4 - Hinged end / Hinged end\n');
fprintf(lsm, '%%BEAM.END.LIBERATION\n');
fprintf(lsm, '4\n');
fprintf(lsm, '1     1  1  1  1  1  1  1  1  1  1  1  1\n');
fprintf(lsm, '2     1  1  1  1  1  0  1  1  1  1  1  1\n');
fprintf(lsm, '3     1  1  1  1  1  1  1  1  1  1  1  0\n');
fprintf(lsm, '4     1  1  1  1  1  0  1  1  1  1  1  0\n');
fprintf(lsm, '\n');

%--------------------------------------------------------------------------
% Elements (Navier or Timoshenko)
if drv.nel > 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Provide elements\n');
    fprintf(lsm, 'First line: Total number of elements\n');
    fprintf(lsm, 'Following lines: Element ID, Material ID, Cross-section ID, End lib. flag,\n');
    fprintf(lsm, '                 Init. Node ID, Final Node ID, vz_X, vz_Y, vz_Z\n');
    if drv.elems(1).type == MEMBER_NAVIER
        fprintf(lsm, '%%ELEMENT.BEAM.NAVIER\n');
    elseif drv.elems(1).type == MEMBER_TIMOSHENKO
        fprintf(lsm, '%%ELEMENT.BEAM.TIMOSHENKO\n');
    end
    
    fprintf(lsm, '%d\n', drv.nel);
    
    for e = 1:drv.nel
        if (drv.elems(e).hingei == CONTINUOUS_END) && (drv.elems(e).hingef == CONTINUOUS_END)
            lib = 1;
        elseif (drv.elems(e).hingei == HINGED_END) && (drv.elems(e).hingef == CONTINUOUS_END)
            lib = 2;
        elseif (drv.elems(e).hingei == CONTINUOUS_END) && (drv.elems(e).hingef == HINGED_END)
            lib = 3;
        elseif (drv.elems(e).hingei == HINGED_END) && (drv.elems(e).hingef == HINGED_END)
            lib = 4;
        end
        
        fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d  %d  %d\n', e,...
                drv.elems(e).material.id, drv.elems(e).section.id, lib,...
                drv.elems(e).nodes(1).id, drv.elems(e).nodes(2).id,...
                drv.elems(e).vz(1), drv.elems(e).vz(2), drv.elems(e).vz(3));
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Nodal prescribed displacements
nnpd = 0;
for n = 1:drv.nnp
    if ~isempty(drv.nodes(n).prescDispl)
        nnpd = nnpd + 1;
    end
end

if nnpd ~= 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Specify nodal prescribed displacements\n');
    fprintf(lsm, 'First line: Total number of nodes with prescribed displacement\n');
    fprintf(lsm, 'Following lines: Node ID, dX [mm], dY [mm], dZ [mm], rX [rad], rY [rad], rZ [rad]\n');
    fprintf(lsm, '%%LOAD.CASE.NODAL.DISPLACEMENT\n');
    fprintf(lsm, '%d\n', nnpd);
    for n = 1:drv.nnp
        if ~isempty(drv.nodes(n).prescDispl)
            fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d\n', n,...
                    1e3*drv.nodes(n).prescDispl(1), 1e3*drv.nodes(n).prescDispl(2),...
                    1e3*drv.nodes(n).prescDispl(3),     drv.nodes(n).prescDispl(4),...
                        drv.nodes(n).prescDispl(5),     drv.nodes(n).prescDispl(6));
        end
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Nodal loads
nnnl = 0;
for n = 1:drv.nnp
    if ~isempty(drv.nodes(n).nodalLoad)
        nnnl = nnnl + 1;
    end
end

if nnnl ~= 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Specify applied nodal loads\n');
    fprintf(lsm, 'First line: Total number of nodes with applied load\n');
    fprintf(lsm, 'Following lines: Node ID, fX [kN], fY [kN], fZ [kN], mX [kNm], mY [kNm], mZ [kNm]\n');
    fprintf(lsm, '%%LOAD.CASE.NODAL.FORCE\n');
    fprintf(lsm, '%d\n', nnnl);
    for n = 1:drv.nnp
        if ~isempty(drv.nodes(n).nodalLoad)
            fprintf(lsm, '%d     %d  %d  %d  %d  %d  %d\n', n,...
                    drv.nodes(n).nodalLoad(1), drv.nodes(n).nodalLoad(2),...
                    drv.nodes(n).nodalLoad(3), drv.nodes(n).nodalLoad(4),...
                    drv.nodes(n).nodalLoad(5), drv.nodes(n).nodalLoad(6));
        end
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Element uniformly distributed loads
neul = 0;
for e = 1:drv.nel
    if ~isempty(drv.elems(e).load.uniformGbl)
        neul = neul + 1;
    end
end

if neul ~= 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Specify element uniformly distributed loads\n');
    fprintf(lsm, 'First line: Total number of elements with uniformly distributed load\n');
    fprintf(lsm, 'Following pair of lines:\n');
    fprintf(lsm, '                         Element ID, Load direction (0->Global, 1->Local)\n');
    fprintf(lsm, '                         Qx [kN/m], Qy [kN/m], Qz [kN/m]\n');
    fprintf(lsm, '%%LOAD.CASE.BEAM.UNIFORM\n');
    fprintf(lsm, '%d\n', neul);
    for e = 1:drv.nel
        if ~isempty(drv.elems(e).load.uniformGbl)
            fprintf(lsm, '%d     %d\n', e, drv.elems(e).load.uniformDir);
            
            if drv.elems(e).load.uniformDir == GLOBAL_LOAD
                fprintf(lsm, '%d  %d  %d\n',...
                        drv.elems(e).load.uniformGbl(1),...
                        drv.elems(e).load.uniformGbl(2),...
                        drv.elems(e).load.uniformGbl(3));
                
            elseif drv.elems(e).load.uniformDir == LOCAL_LOAD
                fprintf(lsm, '%d  %d  %d\n',...
                        drv.elems(e).load.uniformLcl(1),...
                        drv.elems(e).load.uniformLcl(2),...
                        drv.elems(e).load.uniformLcl(3));
            end
        end
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Element linearly distributed loads
nell = 0;
for e = 1:drv.nel
    if ~isempty(drv.elems(e).load.linearGbl)
        nell = nell + 1;
    end
end

if nell ~= 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Specify element linearly distributed loads\n');
    fprintf(lsm, 'First line: Total number of elements with linearly distributed load\n');
    fprintf(lsm, 'Following triple of lines:\n');
    fprintf(lsm, '                         Element ID, Load direction (0->Global, 1->Local)\n');
    fprintf(lsm, '                         Qx_init [kN/m], Qy_init [kN/m], Qz_init [kN/m]\n');
    fprintf(lsm, '                         Qx_final [kN/m], Qy_final [kN/m], Qz_final [kN/m]\n');
    fprintf(lsm, '%%LOAD.CASE.BEAM.LINEAR\n');
    fprintf(lsm, '%d\n', nell);
    for e = 1:drv.nel
        if ~isempty(drv.elems(e).load.linearGbl)
            fprintf(lsm, '%d     %d\n', e, drv.elems(e).load.linearDir);
            
            if drv.elems(e).load.linearDir == GLOBAL_LOAD
                fprintf(lsm, '%d  %d  %d\n',...
                        drv.elems(e).load.linearGbl(1),...
                        drv.elems(e).load.linearGbl(2),...
                        drv.elems(e).load.linearGbl(3));
                fprintf(lsm, '%d  %d  %d\n',...
                        drv.elems(e).load.linearGbl(4),...
                        drv.elems(e).load.linearGbl(5),...
                        drv.elems(e).load.linearGbl(6));
                
            elseif drv.elems(e).load.linearDir == LOCAL_LOAD
                fprintf(lsm, '%d  %d  %d\n',...
                        drv.elems(e).load.linearLcl(1),...
                        drv.elems(e).load.linearLcl(2),...
                        drv.elems(e).load.linearLcl(3));
                fprintf(lsm, '%d  %d  %d\n',...
                        drv.elems(e).load.linearLcl(4),...
                        drv.elems(e).load.linearLcl(5),...
                        drv.elems(e).load.linearLcl(6));
            end
        end
    end
    fprintf(lsm, '\n');
end

%--------------------------------------------------------------------------
% Element temperature variations
netv = 0;
for e = 1:drv.nel
    dtx = drv.elems(e).load.tempVar_X;
    dty = drv.elems(e).load.tempVar_Y;
    dtz = drv.elems(e).load.tempVar_Z;
    if (dtx ~= 0) || (dty ~= 0) || (dtz ~= 0)
        netv = netv + 1;
    end
end

if netv ~= 0
    fprintf(lsm, '-----------------------------------------------------------------------\n');
    fprintf(lsm, 'Specify element thermal loads\n');
    fprintf(lsm, 'First line: Total number of elements with thermal load\n');
    fprintf(lsm, 'Following lines: Element ID, dt_x [oC], dt_y [oC], dt_z [oC]\n');
    fprintf(lsm, 'dt_x       --> element temperature variation on its center of gravity axis\n');
    fprintf(lsm, 'dt_y, dt_z --> element temperature gradient relative to local axes Y and Z\n');
    fprintf(lsm, '               (bottomFaceTempVar - topFaceTempVar)\n');
    fprintf(lsm, '%%LOAD.CASE.BEAM.TEMPERATURE\n');
    fprintf(lsm, '%d\n', netv);
    for e = 1:drv.nel
        dtx = drv.elems(e).load.tempVar_X;
        dty = drv.elems(e).load.tempVar_Y;
        dtz = drv.elems(e).load.tempVar_Z;
        if (dtx ~= 0) || (dty ~= 0) || (dtz ~= 0)
            fprintf(lsm, '%d     %d  %d  %d\n', e, dtx, dty, dtz);
        end
    end
    fprintf(lsm, '\n');
end

fprintf(lsm, '%%END');
end