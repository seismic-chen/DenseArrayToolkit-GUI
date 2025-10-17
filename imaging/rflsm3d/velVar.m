function [du_p,du_s,vavg_p,vavg_s] = velVar(vp,vs,nx,ny,nz)


    du_p = zeros(nz,nx,ny);  
    vavg_p = zeros(nz,1);
    du_s = zeros(nz,nx,ny);
    vavg_s = zeros(nz,1);
    
    for iz = 1:nz
        vp_iz = squeeze(vp(iz,:,:));
        vavg_p(iz) = mean(vp_iz(:)); 
        vs_iz = squeeze(vs(iz,:,:));
        vavg_s(iz) = mean(vs_iz(:)); 
        for ix = 1:nx
            for iy = 1:ny
                du_p(iz,ix,iy) = (1/vp(iz,ix,iy)-1/vavg_p(iz)) + rand(1)*10^-10;            % 扰动速度 纵波
                du_s(iz,ix,iy) = (1/vs(iz,ix,iy)-1/vavg_s(iz)) + rand(1)*10^-10;            % 扰动速度 横波
            end
        end
    end

end