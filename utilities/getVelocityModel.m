function velocityModel = getVelocityModel()
xpad = 50;
xmax = 350;
model = obtain_crust1_QB();
dx = 4;
x = 0-xpad:dx:xmax+xpad;
nx = length(x);
dz = 1;
z = 0:dz:100;
nz = length(z);
vp = interp1(model(:,1),model(:,2),z,'nearest','extrap');
vs = interp1(model(:,1),model(:,3),z,'nearest','extrap');
vel = repmat(vp(:),1,nx);
vel_s = repmat(vs(:),1,nx);

% smooth the velocity model
N=5;
[vel,~]=moving_avg(vel,N,'constant',2);
[vel,~]=moving_avg(vel,N,'constant');
[vel_s,~]=moving_avg(vel_s,N,'constant',2);
[vel_s,~]=moving_avg(vel_s,N,'constant');

velocityModel.x = x;
velocityModel.z = z;
velocityModel.xmax = xmax;
velocityModel.nx = nx;
velocityModel.nz = nz;
velocityModel.dx = dx;
velocityModel.dz = dz;
velocityModel.xpad = xpad;
velocityModel.vp = vel;
velocityModel.vs = vel_s;

figure;
imagesc(x,z,vel); % hold on; scatter(rx,zeros(size(rx)),100,'r^','filled')
end