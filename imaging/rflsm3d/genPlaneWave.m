function dsrc = genPlaneWave(src_func,take_off,back_az,xo,yo,x,y,vp,nt,dt,src_type,fpeak)

v_mean = mean(vp(end,:,:),'all');

theta = deg2rad(take_off);

tbaz = back_az + 180;
phi = deg2rad(tbaz);

nx = length(x);
ny = length(y);

ncoordproj = length(xo);

td = zeros(ncoordproj,1);

coord = zeros(ncoordproj,2);

n1 = sin(theta)*sin(phi);
n2 = sin(theta)*cos(phi);
n3 = cos(theta);

coord(:,1) = xo;
coord(:,2) = yo;

for k = 1:size(coord,1)
    xi = coord(k,1);
    yi = coord(k,2);
    sta = [xi,yi,0];
    kv = [n1,n2,n3];

    dist = (sta*kv')/dot([n1,n2,n3],[n1,n2,n3]);
    dtime = dist / v_mean;
    td(k,1) = dtime;
end

if min(td(:))<0
    td = td + abs(min(td(:)));
end

tmatrix = reshape(td,ny,nx);
timediff = (tmatrix');
%%
dsrc = zeros(nt,nx,ny);
t = (0:nt-1)*dt;
for indx = 1:nx
    for indy = 1:ny
        if strcmp(src_type,'ricker')
            par = pi*fpeak*(t-5);
            src_func = 10*exp(-par.*par).*(1-2*par.*par);  % ricker wavelet
            dsrc(:,indx,indy) = fftShift(src_func',t,timediff(indx,indy));

        else
            dsrc(:,indx,indy) = fftShift(src_func,t,timediff(indx,indy));
        end
    end
end


% ------------------------------------------------------------------------
% % TO CHECK THE TIME
% tmp = rot90(timediff);
% figure
% set(gcf,"Position",[50 50 1300 580],"Color",'white')
% subplot(121)
% imagesc(x,flip(y(:)),tmp)
% axis equal
% xlim([min(x)-30 max(x)+30 ])
% ylim([min(y)-30 max(y)+30 ])
% xlabel('---> X distance')
% ylabel('---> Y distance')
% set(gca,"YDir",'normal');colorbar
% title(['Time Diff: ','Baz ',num2str(back_az),' deg'])
% subplot(122)
% scatter(xo(:),yo(:),80,td(:),'filled','o','MarkerEdgeColor','k')
% axis equal
% xlabel('---> X distance')
% ylabel('---> Y distance')
% set(gca,"YDir",'normal');colorbar
% title(['Time Diff: ','Baz ',num2str(back_az),' deg'])
% xlim([min(xo)-30 max(xo)+30 ])
% ylim([min(yo)-30 max(yo)+30 ])
end
