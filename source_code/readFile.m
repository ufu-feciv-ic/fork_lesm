%% Read file Function
%
% This is an auxiliary file of the <main.html LESM> (Linear Elements
% Structure Model) program that contains functions to read a neutral-format
% file with the _.lsm_ extension 
% This neutral-format file contains all information about a linear element
% structural model.
%
%% Main function
% Output:
%  vs:      flag for version compatibility
%  print:   object of the Print class
% Input arguments:
%  fid: integer identifier of the input file
%  drv: handle to an object of the Drv class
function [vs,print] = readFile(fid,drv)
    % Reads input file
    flag = 1;
    while flag
        % Get file line
        tline = fgetl(fid);
        % Get rid of blank spaces
        string = deblank(tline);
        
        % Look for for tag strings
        switch string
            case '%HEADER.VERSION'
                vs = readVersion(fid);
            case '%HEADER.ANALYSIS'
                [print] = readHeaderAnalysis(fid,drv);
            case '%NODE.COORD'
                readNodeCoord(fid,drv);
            case '%NODE.SUPPORT'
                readNodeSupport(fid,drv);
            case '%MATERIAL.ISOTROPIC'
                readMaterialIsotropic(fid,drv);
            case '%SECTION.PROPERTY'
                readSectionProperty(fid,drv);
            case '%BEAM.END.LIBERATION'
                rotlib = readBeamEndLiberation(fid);
            case '%ELEMENT.BEAM.NAVIER'
                readElementBeamNavier(fid,drv,rotlib);
            case '%ELEMENT.BEAM.TIMOSHENKO'
                readElementBeamTimoshenko(fid,drv,rotlib);
            case '%LOAD.CASE.NODAL.FORCE'
                readNodeAppliedForceMoment(fid,drv);
                case '%LOAD.CASE.NODAL.DISPLACEMENT'
                readNodePrescDisplRot(fid,drv);
            case '%LOAD.CASE.BEAM.UNIFORM'
                readElemAppliedUnifLoad(fid,drv);
            case '%LOAD.CASE.BEAM.LINEAR'
                readElemAppliedLinearLoad(fid,drv);
            case '%LOAD.CASE.BEAM.TEMPERATURE'
                readElemThermalLoad(fid,drv);
        end
        
        if strcmp(string,'%END')
            flag = 0;
        end
        
    end
    
    % Close file
    fclose(fid);
end

%% Auxiliary functions
%--------------------------------------------------------------------------
function vs = readVersion(fid)
    % Read program version number
    version = fscanf(fid,'%f',1);
    
    % Check version compatibility
    if (version == 1.0) || (version == 1.1)
        vs = 1;
    else
        vs = 0;
        return
    end
end

%--------------------------------------------------------------------------
function [print] = readHeaderAnalysis(fid,drv)
    % Get file line
    tline = fgetl(fid);
    % Get rid of blank spaces
    string = deblank(tline);
    
    % Get target analysis type
    switch string
        case '''TRUSS2D'''
            drv.anm = Anm_Truss2D();
            print = Print_Truss2D(drv);
        case '''FRAME2D'''
            drv.anm = Anm_Frame2D();
            print = Print_Frame2D(drv);
        case '''GRILLAGE'''
            drv.anm = Anm_Grillage();
            print = Print_Grillage(drv);
        case '''TRUSS3D'''
            drv.anm = Anm_Truss3D();
            print = Print_Truss3D(drv);
        case '''FRAME3D'''
            drv.anm = Anm_Frame3D();
            print = Print_Frame3D(drv);
    end
end

%--------------------------------------------------------------------------
function readNodeCoord(fid,drv)
    % Read and store total number of nodes and equations
    nnp = fscanf(fid,'%f',1);
    drv.nnp = nnp;
    drv.neq = nnp * drv.anm.ndof;
    
    % Initialize vector of objects of the Node class
    nodes(1,nnp) = Node();
    drv.nodes = nodes;
    
    for n = 1:nnp
        b = fscanf(fid,'%f',4);
        % b(1) = node id
        % b(2) = X coordinate
        % b(3) = Y coordinate
        % b(4) = Z coordinate
        drv.nodes(n).id = b(1);
        drv.nodes(n).coord = [b(2) b(3) b(4)];
    end
end

%--------------------------------------------------------------------------
function readNodeSupport(fid,drv)
    % Read total number of nodes with essential boundary conditions (supports)
    nnebc = fscanf(fid,'%f',1);
    
    for n = 1:nnebc
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Dx
        % b(3) = Dy
        % b(4) = Dz
        % b(5) = Rx
        % b(6) = Ry
        % b(7) = Rz
        drv.nodes(b(1)).ebc = [b(2) b(3) b(4) b(5) b(6) b(7)];
    end
end

%--------------------------------------------------------------------------
function readMaterialIsotropic(fid,drv)
    % Read and store total number of materials
    nmat = fscanf(fid,'%f',1);
    drv.nmat = nmat;
    
    % Initialize vector of objects of the Material class
    materials(1,nmat) = Material();
    drv.materials = materials;
    
    for m = 1:nmat
        b = fscanf(fid,'%f',4);
        % b(1) = material id
        % b(2) = elasticity modulus [MPa]
        % b(3) = possion ratio
        % b(4) = thermal expansion coefficient  [/oC]
        drv.materials(m) = Material(b(1),1e3*b(2),b(3),b(4));
    end
end

%--------------------------------------------------------------------------
function readSectionProperty(fid,drv)
    % Read and store total number of cross-sections
    nsec = fscanf(fid,'%f',1);
    drv.nsec = nsec;
    
    % Initialize vector of objects of the Section class
    sections(1,nsec) = Section();
    drv.sections = sections;
    
    for s = 1:nsec
        b = fscanf(fid,'%f',[1 9]);
        % b(1) = Cross-section id
        % b(2) = Ax [cm2]
        % b(3) = Ay [cm2]
        % b(4) = Az [cm2]
        % b(5) = Ix [cm4]
        % b(6) = Iy [cm4]
        % b(7) = Iz [cm4]
        % b(8) = Hy [cm]
        % b(9) = Hz [cm]
        drv.sections(s).id = b(1);
        drv.sections(s).area_x = 1e-4*b(2);
        drv.sections(s).area_y = 1e-4*b(3);
        drv.sections(s).area_z = 1e-4*b(4);
        drv.sections(s).inertia_x = 1e-8*b(5);
        drv.sections(s).inertia_y = 1e-8*b(6);
        drv.sections(s).inertia_z = 1e-8*b(7);
        drv.sections(s).height_y = 1e-2*b(8);
        drv.sections(s).height_z = 1e-2*b(9);
    end
end

%--------------------------------------------------------------------------
function rotlib = readBeamEndLiberation(fid)
    include_constants;
    
    % Read number of beam end rotation liberation properties
    nrotlib = fscanf(fid,'%f',1);
    rotlib = zeros(2,nrotlib);
    
    % Read beam end rotation liberation (hinge) data
    for i = 1:nrotlib
        b = fscanf(fid,'%f',[1 13]);
        % b( 1) = rotation liberation property id (not used)
        % b( 2) = idx
        % b( 3) = idy
        % b( 4) = idz
        % b( 5) = irx
        % b( 6) = iry
        % b( 7) = irz
        % b( 8) = jdx
        % b( 9) = jdy
        % b(10) = jdz
        % b(11) = jrx
        % b(12) = jry
        % b(13) = jrz
        % If rotation liberation in any of the directions is free, it is
        % considered that it is free in all directions, thus it will receive a
        % hinged flag
        if (b(5)==0) || (b(6)==0) || (b(7)==0)
            rotlib(1,i) = HINGED_END;
        else
            rotlib(1,i) = CONTINUOUS_END;
        end
        
        if (b(11)==0) || (b(12)==0) || (b(13)==0)
            rotlib(2,i) = HINGED_END;
        else
            rotlib(2,i) = CONTINUOUS_END;
        end
    end
end

%--------------------------------------------------------------------------
function readElementBeamNavier(fid,drv,rotlib)
    % Read and store total number of elements
    nel = fscanf(fid,'%f',1);
    drv.nel = nel;
    
    % Initialize vector of objects of the Elem_Navier sub-class
    elems(1,nel) = Elem_Navier();
    drv.elems = elems;
    
    for e = 1:nel
        b = fscanf(fid,'%f',[1 9]);
        % b(1) = element id
        % b(2) = material id
        % b(3) = section property id
        % b(4) = end release id
        % b(5) = init node id
        % b(6) = final node id
        % b(7) = vz_X
        % b(8) = vz_Y
        % b(9) = vz_Z
        drv.elems(e) = Elem_Navier(drv.anm,...
                                   drv.materials(b(2)),...
                                   drv.sections(b(3)),...
                                   [drv.nodes(b(5)), drv.nodes(b(6))],...
                                   rotlib(1,b(4)), rotlib(2,b(4)),...
                                   [b(7), b(8), b(9)],...
                                   [],[],[],[],[],[],[],[]);
    end
end

%--------------------------------------------------------------------------
function readElementBeamTimoshenko(fid,drv,rotlib)
    % Read and store total number of elements
    nel = fscanf(fid,'%f',1);
    drv.nel = nel;
    
    % Initialize vector of objects of the Elem_Timoshenko class
    elems(1,nel) = Elem_Timoshenko();
    drv.elems = elems;
    
    for e = 1:nel
        b = fscanf(fid,'%f',[1 9]);
        % b(1) = element id
        % b(2) = material id
        % b(3) = section property id
        % b(4) = end release id
        % b(5) = init node id
        % b(6) = final node id
        % b(7) = vz_X
        % b(8) = vz_Y
        % b(9) = vz_Z
        
        % Create element object
        drv.elems(e) = Elem_Timoshenko(drv.anm,...
                                       drv.materials(b(2)),...
                                       drv.sections(b(3)),...
                                       [drv.nodes(b(5)), drv.nodes(b(6))],...
                                       rotlib(1,b(4)), rotlib(2,b(4)),...
                                       [b(7), b(8), b(9)],...
                                       [],[],[],[],[],[],[],[]);
    end
end

%--------------------------------------------------------------------------
function readNodeAppliedForceMoment(fid,drv)
    % Read total number of nodal loads
    n_nodalload = fscanf(fid,'%f',1);
    
    for i = 1:n_nodalload
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Fx [kN]
        % b(3) = Fy [kN]
        % b(4) = Fz [kN]
        % b(5) = Mx [kNm]
        % b(6) = My [kNm]
        % b(7) = Mz [kNm]
        drv.nodes(b(1)).nodalLoad = [b(2) b(3) b(4) b(5) b(6) b(7)];
    end
end

%--------------------------------------------------------------------------
function readNodePrescDisplRot(fid,drv)
    % Read number of nodes with prescribed displacements
    nnprescdispl = fscanf(fid,'%f',1);
    
    for i = 1:nnprescdispl
        b = fscanf(fid,'%f',[1 7]);
        % b(1) = node id
        % b(2) = Dx [mm]
        % b(3) = Dy [mm]
        % b(4) = Dz [mm]
        % b(5) = Rx [rad]
        % b(6) = Ry [rad]
        % b(7) = Rz [rad]
        drv.nodes(b(1)).prescDispl(1) = 1e-3*b(2);
        drv.nodes(b(1)).prescDispl(2) = 1e-3*b(3);
        drv.nodes(b(1)).prescDispl(3) = 1e-3*b(4);
        drv.nodes(b(1)).prescDispl(4) = b(5);
        drv.nodes(b(1)).prescDispl(5) = b(6);
        drv.nodes(b(1)).prescDispl(6) = b(7);
    end
end

%--------------------------------------------------------------------------
function readElemAppliedUnifLoad(fid,drv)
    % Read number of elements with applied uniformely distrib. load
    n_uniformload = fscanf(fid,'%f',1);
    
    for i = 1:n_uniformload
        b = fscanf(fid,'%f',[1 2]);
        % b(1) = element id
        % b(2) = direction (global or local system)
        drv.elems(b(1)).load.uniformDir = b(2);
        
        c = fscanf(fid,'%f',[1 3]);
        % c(1) = Qx [kN/m]
        % c(2) = Qy [kN/m]
        % c(3) = Qz [kN/m]
        drv.elems(b(1)).load.setUnifLoad([c(1),c(2),c(3)],b(2));
    end
end

%--------------------------------------------------------------------------
function readElemAppliedLinearLoad(fid,drv)
    % Read number of elements with applied linearly distrib. load
    n_linearload = fscanf(fid,'%f',1);
    
    for i = 1:n_linearload
        b = fscanf(fid,'%f',[1 2]);
        % b(1) = element id
        % b(2) = direction (local or global system)
        drv.elems(b(1)).load.linearDir = b(2);
        
        c = fscanf(fid,'%f',[1 3]);
        % c(1) = Qxi [kN/m]
        % c(2) = Qyi [kN/m]
        % c(3) = Qzi [kN/m]
        d = fscanf(fid,'%f',[1 3]);
        % d(1) = Qxj [kN/m]
        % d(2) = Qyj [kN/m]
        % d(3) = Qzj [kN/m]
        drv.elems(b(1)).load.setLinearLoad([c(1),c(2),c(3),d(1),d(2),d(3)],b(2));
    end
end

%--------------------------------------------------------------------------
function readElemThermalLoad(fid,drv)
    % Read number of elements with thermal load
    n_tempvar = fscanf(fid,'%f',1);
    
    for i = 1:n_tempvar
        b = fscanf(fid,'%f',[1 4]);
        % b(1) = element id
        % b(2) = dtx [oC]
        % b(3) = dty [oC]
        % b(4) = dtz [oC]
        drv.elems(b(1)).load.tempVar_X = b(2);
        drv.elems(b(1)).load.tempVar_Y = b(3);
        drv.elems(b(1)).load.tempVar_Z = b(4);
    end
end