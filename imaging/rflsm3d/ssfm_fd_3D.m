function [mod,mod_source,mod_receiver] = ssfm_fd_3D(img,planewave,save_wavefield,param,if_cg)

    % parameters
    
    vp = param.vp;
    vs = param.vs;
    nt = param.nt;
    dt = param.dt;
    dx = param.dx;
    dy = param.dy;
    dz = param.dz;
    flow = param.flow;
    fhigh = param.fhigh;
    bc = param.bc;
    
    %-------------------------------------------------------------------------
    %-------------------------------------------------------------------------
    if ~exist('if_cg','var')
      if_cg = 0;
    end
    
    [nz,nx,ny] = size(vp);
    
    if if_cg==1
       img = reshape(img,nz,nx,ny);
    end
    
    nf  = 2^nextpow2(nt);
    nsx = 2^nextpow2(nx);
    nsy = 2^nextpow2(ny);
    
    w = 2*pi/(nf*dt)*[0:(nf/2-1),-(nf/2):-1]; 
    
    ind = (w==0);
    w(ind)=1E-30;
    
    iw1 = floor(flow*dt*nf)+1;
    iw2 = floor(fhigh*dt*nf)+1;
    % if frequency exceeding the nyquest freq
    if iw2 > floor(nf/2)+1
        iw2 = floor(nf/2)+1;
    end
    
    %% vel pertu
    [du_p,du_s,vavg_p,vavg_s] = velVar(vp,vs,nx,ny,nz);
    
    %% plane source
    dsc = fft(planewave,nf,1);
   
    % output data in frequence domian
    outf = zeros(nf,nsx,nsy);
    outf_source = zeros(nz,nf,nsx,nsy);
    outf_receiver = zeros(nz,nf,nsx,nsy);
    
    parfor iw = iw1:iw2
        source_input = reshape(dsc(iw,:,:),nx,ny);
        % [swave,swave_full] = sspropag_op_3D(source_input,vavg_p,du_p,nx,dx,nz,dz,ny,dy,w(iw),1,bc,'source',save_wavefield);
        [swave,swave_full] = sspropag_op_3D_fast(source_input,vavg_p,du_p,nx,dx,nz,dz,ny,dy,w(iw),1,bc,'source',save_wavefield);
        % swave(nz,nx,ny) F-X domain for each frequency
        % swave_full(nz,nsx,nsy) F-K domain for each frequency
        
        receiver_input = (img.*swave);
        % [rwave_s,rwave_full_s] = sspropag_op_3D(receiver_input,vavg_s,du_s,nx,dx,nz,dz,ny,dy,w(iw),1,bc,'receiver',save_wavefield);
        [rwave_s,rwave_full_s] = sspropag_op_3D_fast(receiver_input,vavg_s,du_s,nx,dx,nz,dz,ny,dy,w(iw),1,bc,'receiver',save_wavefield);
        % rwave_s(nsx,nsy) : f-k domain of each iw ,iz = 1;
        % rwave_full_s(nz,nsx,nsy): f-k domain of each iw
        outf(iw,:,:) = rwave_s;

        
        if save_wavefield
            outf_source(:,iw,:,:) = swave_full;
            % outf_receiver(:,iw,:,:) = rwave_full_s;        % whole field           
            % outf_receiver(:,iw,:)=rwave_full_s+rwave_full_p; % only used
            % when simulating P phase in RF
        end  
    end

    %--------------transform outf to tx domain-------------
    % iz = 1;
    r = real(ifft(ifft2(outf,nf,nsx),nsy,3));
    modtemp = r(1:nt,1:nx,1:ny);  

    mod = modtemp;
    
    if if_cg == 1
        mod = modtemp(:);
    end
    
    if save_wavefield 
        % surface
        tmp_source = reshape(outf_source(1,:,:,:),nf,nsx,nsy);
        % tmp_receiver = reshape(outf_receiver(1,:,:,:),nf,nsx,nsy);
        s_full = real(ifft(ifft2(tmp_source,nf,nsx),nsy,3));
        % r_full = real(ifft(ifft2(tmp_receiver,nf,nsx),nsy,3));
        mod_source = s_full(1:nt,1:nx,1:ny);
        % mod_receiver(iz,:,:,:) = r_full(1:nt,1:nx,1:ny);
    else
        mod_source=[];
    end
    mod_receiver=[];

end
