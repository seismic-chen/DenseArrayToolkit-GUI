function [mig,dp] = runMigration(rfshift,take_off,back_az,src_func,save_wavefield,gridStruct,param)

    param = param.paramMig;

    xo = gridStruct.XInOriginalCoord(:); 
    yo = gridStruct.YInOriginalCoord(:);

    x = gridStruct.x;             
    y = gridStruct.y;

    dt = param.dt;
    nt = param.nt;
    vp = param.vp;
    src_type = param.src_type;
    fpeak = param.fpeak;

    ispred = param.ispred;
    
    %% 
    dsrc = genPlaneWave(src_func,take_off,back_az,xo,yo,x,y,vp,nt,dt,src_type,fpeak);
    
    disp('==========>')
    disp('do migration')
    mig = ssfm_adj_3D(rfshift,dsrc,save_wavefield,param);

    % predict data
    if ispred
        save_wavefield = 0;  % 1: only to calculate time diff
        [dp,~,~] = ssfm_fd_3D(mig,dsrc,save_wavefield,param);
        S = repmat(any(rfshift),size(rfshift,1),1);

        dp = S.*dp;
    else
        dp = [];
    end


end