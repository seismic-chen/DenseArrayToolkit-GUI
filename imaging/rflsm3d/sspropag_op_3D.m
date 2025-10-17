function [wavf,wavf_full] = sspropag_op_3D(in,vavg,du,nx,dx,nz,dz,ny,dy,w,iflag,bc,wavf_type,save_wavefield)
% ssf 主要程序
% 包含震源波场和检波点波场计算
% 2023.3.30

% phase shift wave propogator with split-step correction for 1 frequency
% note that it is only used to extrapolate from 1 to nz (not nz to 1)

% if iflag == -1 do migration propogation
% if iflag == 1  do modeling propogation

% if dz > 0 apply downward propogation
% if dz < 0 apply upward propogation

% input:  in   -- data/model to be propogated
%          w   -- angular frequency
%       vavg   -- average velocity of each depth
%         du   -- purtubation of slowness (size of vel model)
%         dx   -- x spatial interval
%         dz   -- z spatial interval
%         dt   -- time interval
%  wavf_type   -- wavefiled type (source or receiver)
%  save_wavefield -- set to 1 to save wavefields at each depth
% output:
%         wavf -- propogated wavefield at a given frequency
%                with split step corrections
%      wav_full -- 
%% 边界条件 
if bc == 1
    if nx == 1
        Lx = 0;
    else
        Lx = 5;    
    end
    ix = 1:1:Lx;
    tap1 = exp(-(0.1*(Lx-ix)).^2);
    wex = [tap1,ones(1,nx-2*Lx),fliplr(tap1)];
    if ny == 1
        Ly = 0;
    else
        Ly = 5;
    end
    iy = 1:1:Ly;
    tap2 = exp(-(0.1*(Ly-iy)).^2);
    wey = [tap2,ones(1,ny-2*Ly),fliplr(tap2)];

else
    wex = ones(1,nx);
    wey = ones(1,ny);
end
wex = repmat(wex,ny,1);
wey = repmat(wey,nx,1);
we = (wex'+wey)/2;

%% ---------------------------------------------------------------%
nsx = 2^nextpow2(nx);
nsy = 2^nextpow2(ny);
if nsx == 1
    Kx = 0;
else
    Kx = 2*pi/(dx*nsx)*[0:(nsx/2-1),-(nsx/2):-1]; 
end
if nsy == 1
    Ky = 0;
else
    Ky = 2*pi/(dy*nsy)*[0:(nsy/2-1),-(nsy/2):-1];
end
%% ---------------------------------------------------------------%
if save_wavefield
    wavf_full = zeros(nz,nsx,nsy);         % 
else
    wavf_full = [];
end
% 
if iflag == -1 
    % in: fx domain
    M0_old = fft2(in,nsx,nsy);
    wavf = zeros(nz,nx,ny);

    for iz = 1:nz

        S0 = zeros(nsx,nsy);
        for ikx = 1:nsx
            for iky = 1:nsy

                arg = w^2/vavg(iz)^2 - Kx(ikx)^2- Ky(iky)^2;    
                
                if arg >= 0
                    S0(ikx,iky) = exp(1i*(sqrt(arg))*dz); 
                else
                    S0(ikx,iky) = exp(-1*(sqrt(-arg))*abs(dz));
                end
            end
        end

        M0 = M0_old.*S0;
        P0 = ifft2(M0,nsx,nsy);

        C0 = exp(1i*w*dz * reshape(du(iz,:,:),nx,ny));

        P1 = P0(1:nx,1:ny).*C0;
 
        P1 = P1.* we;        % z+1 

        wavf(iz,:,:) = P1; 

        M0_old = fft2(P1,nsx,nsy); 
        
        if save_wavefield
            % surface
            % if iz == 1
            wavf_full(iz,:,:) = M0_old;
            % end
        end        
    end
end

%  forward modeling
if iflag == 1
    switch wavf_type
        case 'source'
            % dshift in fx domain ---> fk domain
            mfk = fft2(in,nsx,nsy);      % (x,y) --> (kx,ky) for each w at depth zi
            wavf = zeros(nz,nx,ny);         
        case 'receiver'
            % in: (nz,nx,ny) domain for each frequency
            mfk = zeros(nz,nsx,nsy);
            % fx --> fk
            for iz = 1:nz
                mfk(iz,:,:) = fft2(reshape(in(iz,:,:),nx,ny),nsx,nsy);
            end
            wavf = zeros(nsx,nsy);
    end

    for iz = nz:-1:1
        S0 = zeros(nsx,nsy);    % 
        for ikx = 1:nsx
            for iky = 1:nsy
                arg = w^2/vavg(iz)^2 - Kx(ikx)^2 - Ky(iky)^2;    % vertical wavenumber kz
                if arg >= 0
                    S0(ikx,iky) = exp(-1i*(sqrt(arg))*dz);       % phase shift operator
                else
                    S0(ikx,iky) = exp(-1*(sqrt(-arg))*abs(dz));
                end
            end
        end

        %% 
        switch wavf_type
            case 'source'
                if iz == nz
                    mfk_old = mfk;    % (kx,ky) domain, initial condition at the bottom of the model
                else
                    % fk domain * operater
                    mfk = mfk_old.*S0;
                    %  fk ---> fx
                    P0 = ifft2(mfk,nsx,nsy);
                    % fx domain * operater
                    C0 = exp(-1i*w*(dz) * reshape(abs(du(iz,:,:)),nx,ny)); 
                    % C0 = ones(nx,ny);
                    P1 = P0(1:nx,1:ny).*C0;
                    P1 = P1.*we;

                    wavf(iz,:,:) = P1;    % (x,y) domain at depth zi
                    mfk_old = fft2(P1,nsx,nsy);

                end
                if save_wavefield
                    wavf_full(iz,:,:) = mfk_old;
                end

            case 'receiver'
                % fk -- > fx
                P0 = ifft2((wavf + reshape(mfk(iz,:,:),nsx,nsy)),nsx,nsy);
                C0 = exp(-1i*w*dz * reshape(du(iz,:,:),nx,ny));

                P1 = P0(1:nx,1:ny).*C0;     % (nx,ny)
                P1 = P1.* we;
                
                % FX---->FK
                wavf = fft2(P1,nsx,nsy).*S0;
                if save_wavefield
                    wavf_full(iz,:,:) = wavf;
                end
                
        end    % end switch     
    end     % end iz
end     % end iflag
end     % end function
