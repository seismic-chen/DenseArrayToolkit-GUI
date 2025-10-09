function param = setMigParam3D(gridStruct)

    param.rx = gridStruct.rx;
    param.ry = gridStruct.ry;
    
    param.dx = gridStruct.dx;
    param.dy = gridStruct.dy;
    param.dz = gridStruct.dz;
    
    param.x = gridStruct.x;
    param.y = gridStruct.y;
    param.z = gridStruct.z;
    param.nx = gridStruct.nx;
    param.ny = gridStruct.ny;
    param.nz = gridStruct.nz;
    
    param.flow = 0.01;
    param.fhigh = 2.5;
    
    param.src_type = 'P';
    param.fpeak = 1.2;   % do not use, for ricker wave source
    
    param.bc = 1;
    param.plot = 0;

    param.ispred = 0;

    %% velocity model
    param.vp = gridStruct.VP;
    param.vs = gridStruct.VS;

    %% reconstruction or not
    param.isReconRFs = 0;
    
end