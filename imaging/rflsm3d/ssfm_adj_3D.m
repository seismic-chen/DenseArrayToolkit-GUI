function [dmig] = ssfm_adj_3D(din,dsrc,save_wavefield,param,if_cg)

    % parameters
    rx = param.rx;
    ry = param.ry;
    x = param.x;
    y = param.y;
    z = param.z;
    
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
    isplot = param.plot;
      
    %% 
    if ~exist('if_cg','var')
      if_cg = 0;
    end
    
    [nz,nx,ny] = size(vp);
    
    if if_cg==1
       din = reshape(din,nt,nx,ny);
    end
    
    dmig = zeros(nz,nx,ny);
    
    nf = 2^nextpow2(nt);
    nsx = 2^nextpow2(nx);
    nsy = 2^nextpow2(ny);
    w = 2*pi/(nf*dt)*[0:(nf/2-1),-(nf/2):-1];
    ind= (w==0);
    w(ind)=1E-30;
    
    %define loop over frequencies
    iw1 = floor(flow*dt*nf)+1;
    iw2 = floor(fhigh*dt*nf)+1;
    
    if iw2 > floor(nf/2)+1
        iw2=floor(nf/2)+1;
    end
    
    %% vel var
    
    [du_p,du_s,vavg_p,vavg_s] = velVar(vp,vs,nx,ny,nz);
    %%      
    % source
    dsc = fft(dsrc,nf,1);   % time to frequency
    
    % shifted rf
    dfx = fft(din,nf,1);
    
    % image
    img = zeros(nz,nx,ny);
    outf_source = zeros(nz,nf,nsx,nsy);
    outf_receiver = zeros(nz,nf,nsx,nsy);
    parfor iw = iw1:iw2
        source_input = dsc(iw,:,:);
        source_input = reshape(source_input,nx,ny);
        [swave,swave_full] = sspropag_op_3D_fast(source_input,vavg_p,du_p,nx,dx,nz,dz,ny,dy,w(iw), 1,bc,'source', save_wavefield);
        % [swave,swave_full] = sspropag_op_3D(source_input,vavg_p,du_p,nx,dx,nz,dz,ny,dy,w(iw), 1,bc,'source', save_wavefield);
            % swave(nz,nx,ny) F-X domain in frequency slice
        
        % real data
        receiver_input = reshape(dfx(iw,:,:),nx,ny);
        [rwave,rwave_full] = sspropag_op_3D_fast(receiver_input,vavg_s,du_s,nx,dx,nz,dz,ny,dy,w(iw),-1,bc,'receiver',save_wavefield);
        % [rwave,rwave_full] = sspropag_op_3D(receiver_input,vavg_s,du_s,nx,dx,nz,dz,ny,dy,w(iw),-1,bc,'receiver',save_wavefield);
        
        % (nz,nx,ny) 
        %apply imaging condition
        img = img + real(rwave.*conj(swave));
    
        if save_wavefield
            outf_source(:,iw,:,:) = swave_full;
            outf_receiver(:,iw,:,:) = rwave_full; 
            
        end
        
    end
    dmig = img./nf;   
    
    if if_cg==1
        dmig = dmig(:); 
    end
    
    
    if save_wavefield
        mod_source = zeros(nz,nt,nx,ny);
        mod_receiver = zeros(nz,nt,nx,ny); 
        
        parfor iz = 1:nz
            tmp_source = reshape(outf_source(iz,:,:,:),nf,nsx,nsy);
            tmp_receiver = reshape(outf_receiver(iz,:,:,:),nf,nsx,nsy);
    
            s_full = real(ifft(ifft2(tmp_source,nf,nsx),nsy,3));
            r_full = real(ifft(ifft2(tmp_receiver,nf,nsx),nsy,3));
    
            mod_source(iz,:,:,:) = s_full(1:nt,1:nx,1:ny);
            mod_receiver(iz,:,:,:) = r_full(1:nt,1:nx,1:ny);
    
        end
    end
    
    
    if isplot && save_wavefield
        load('roma.mat');
        xline = 150;
        yline = 100;
        rxmin = min(rx);
        rxmax = max(rx);
        rymin = min(ry);
        rymax = max(ry);
        Xl = [xline, xline];
        Yl = [yline, yline];
        Zl = [0, max(z)];
        [X,Y,Z] = meshgrid(x,y,z); 
        % 
        figure
        slice(X,Y,Z,permute(dmig./max(dmig(:)),[3 2 1]),[xline ],[yline],[]);
        set(gca,'ZDir','reverse');
        shading interp;
        axis tight;
        clim([-1 1]);     
        xlim([rxmin rxmax]);
        ylim([rymin rymax]);
        colormap(flipud(roma))
    
        hold on;
        slice(X,Y,Z,permute(param.vp,[3 2 1]),[xline ],[yline],[]);

        pic_num = 1;
        for tt = 81:5:nt-550+1
            swave = reshape(mod_source(:,tt,:,:),nz,nx,ny)./max(mod_source(:));
            rwave = reshape(mod_receiver(:,581 - tt,:,:),nz,nx,ny)./max(mod_receiver(:));
        
            figure(1);
            set(gcf,'position',[102,400,1350,300],'color','white');

            subplot(121)
            h1 = slice(X,Y,Z,permute(swave,[3 2 1]),[xline],[yline],[]);
            title(['Forward (',num2str((tt-1)*dt),'s)']);
            xlabel('x (km)');ylabel('y (km)');zlabel('Depth (km)');
            set(h1,'facealpha',0.8);
            set(gca,'ZDir','reverse');
            shading interp; axis tight; clim([-0.8 0.8]);
            xlim([min(x) max(x)]); ylim([min(y) max(y)]);
            colormap(flipud(roma))
            hold on
            plot3(Xl,Yl,Zl,'LineWidth',1.2,'Color',[0.8 0.8 0.8],'LineStyle','-')
            hold off;
            set(gca,'FontSize',16); 

            subplot(122)
            h2 = slice(X,Y,Z,permute(rwave,[3 2 1]),[xline],[yline],[]);
            title(['Backward (',num2str((tt-1)*dt),'s)']);
            set(h2,'facealpha',0.8);
            set(gca,'ZDir','reverse');
            xlabel('x (km)');ylabel('y (km)');zlabel('Depth (km)');
            shading interp; axis tight; clim([-0.8 0.8]);
            xlim([min(x) max(x)]); ylim([min(y) max(y)]);
            colormap(flipud(roma))
            hold on
            plot3(Xl,Yl,Zl,'LineWidth',1.2,'Color',[0.8 0.8 0.8],'LineStyle','-')
            hold off;
            set(gca,'FontSize',16); 
        
            drawnow;
            F = getframe(gcf);
            I = frame2im(F);
            [I,map] = rgb2ind(I,256);
            if pic_num == 1
                imwrite(I,map,'migration3d.gif','gif','Loopcount',inf,'DelayTime',0.05);
            else
                imwrite(I,map,'migration3d.gif','gif','WriteMode','append','DelayTime',0.05);
            end
            pic_num = pic_num + 1;
        end 
    
    end

return

