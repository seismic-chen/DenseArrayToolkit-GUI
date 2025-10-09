function m=preconT_x_3d(u,nx,ny,nz,dx,dy,dz)


% 
u = reshape(u,nz,nx,ny);
m = zeros(nz,nx,ny);
type = 'pre';

switch type
    case 'reg'
        order = 2;
        if order ==1
            D1 = (diag(ones(1,nx))+diag(-1*ones(1,nx-1),-1)) / dx;
            P = D1^-1;
        elseif order == 2
            D2 = (diag(2*ones(1,nx))+diag(-1*ones(1,nx-1),-1)+diag(-1*ones(1,nx-1),1)) /dx^2;
            P = D2^-1;
        
        end     
        for k = 1:nz
            
            tmp = reshape(u(k,:,:),nx,ny);
            m(k,:,:) =P' * tmp;
        
        end

    case 'pre'
        P=1/2*[1,2,1];
        a1 = zeros(nx,ny);
        a2 = a1;
        for k = 1:nz
            utmp = reshape(u(k,:,:),nx,ny);
            for i = 1:nx
                a1(i,:) = conv(utmp(i,:),P,'same');
            end
            for j = 1:ny
                a2(:,j) = conv(a1(:,j)',P,'same')';
            end
            m(k,:,:) = a2;
        end
end
m = m(:);



