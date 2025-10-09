function [rfshift,src_func,mask] = shiftRFs(rf0,take_off,back_az,xo,yo,x,y,rx,ry,param)

    param = param.paramMig;

    TIME = param.Ti;
    dt = param.dt;
    nt = param.nt;
    nx = param.nx;
    ny = param.ny;
    nz = param.nz;
    vp = param.vp;
    src_type = param.src_type;
    fpeak = param.fpeak;

    dx = param.dx;
    dy = param.dy;
    
    isplot = param.plot;
    
    %% prepeocess and extract wavelet
    [rf1,src_func] = preprocrf(rf0,param);
    rfshift0 = zeros(nt,nx,ny);     % only for plotting
    [rftmp0,mask] = doBinning(rf0, rx, ry, x, y,nx,ny, dx, dy);
    
    %% 
    if param.isReconRFs
        rftmp = param.dotg;
    else
        rfshift = zeros(nt,nx,ny);
        [rftmp,~] = doBinning(rf1, rx, ry, x, y, nx,ny, dx, dy);
    end
    
    %% generate plane wave
    dsrc = genPlaneWave(src_func,take_off,back_az,xo,yo,x,y,vp,nt,dt,src_type,fpeak);
    
    %%  forward
    save_wavefield = 1;
    img = zeros(nz,nx,ny);
    [~,mod_source,~] = ssfm_fd_3D(img,dsrc,save_wavefield,param);
    
    %% apply time diff to rf using cross correlation
    % cross correlate P wave with RF
    it = (0:nt-1)*dt;
    for i = 1:nx
        for j = 1:ny
            sig1 = src_func;
            sig2 = mod_source(:,i,j);
        
%             % cross-correlation
            xc=xcorr(sig2,sig1);
            tax = [-(nt-1):(nt-1)]*dt;
            [~,ind]=max(xc);
            tdelay = tax(ind);
    
            % phase only cross correlation
%             [tdelay, ~, ~] = phase_only_correlation_simple(sig1, sig2, dt);

            % shift receiver function which locate at the regular grid
            rfshift0(:,i,j) = fftShift(rftmp0(:,i,j),it,tdelay);
            rfshift(:,i,j) = fftShift(rftmp(:,i,j),it,tdelay);
        end
    end
    rfshift0(isnan(rfshift0))=0;
    rfshift(isnan(rfshift))=0;
    
    %% plot
    if isplot
        xtrace = 11;
        itrf = it - 10;   % previous to P 10s
        tb = round(5/dt);
        tl = round(50/dt);
        rfplot = rfshift0(:,:,xtrace);
        xprof = rfshift(:,:,xtrace);
        
        figure
        set(gcf,'position',[100 100 1000 400],'color','white')
        subplot(121)
        wigb(rfplot(tb:tl,:),0.8,x,itrf(tb:tl))
        title('shifted RFs')
        xlabel('Distance (km)');
        ylabel('Time (s)');
        set(gca,'fontsize',18);
        box on
    
        subplot(122)
        wigb(xprof(tb:tl,:),0.8,x,itrf(tb:tl))
        title('shifted RFs (tapered)')
        xlabel('Distance (km)');
        ylabel('Time (s)');
        set(gca,'fontsize',18);
        box on
    
    end

end