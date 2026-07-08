clear;clc

LatMin = 40.5;
LatMax = 43;
ynpts = 25;

LonMin = 108;
LonMax = 112;
xnpts = 41;

% Get AK135 reference model
[z, rho, vp, vs, qk, qm] = ak135('cont');
zmax = 100;  % ccp 800km
dz = 1;

% Initialize arrays for storing 3D model data
X = []; Y = []; Z = []; VP = []; VS = [];
knode = 0;

% Create lat/lon grid for interpolation
lonall = linspace(LonMin, LonMax, xnpts);
latall = linspace(LatMin, LatMax, ynpts);

vpgrid = zeros(zmax/dz,xnpts,ynpts);
vsgrid = zeros(zmax/dz,xnpts,ynpts);

% Loop through each grid point to build 3D model
for i = 1:length(latall)
    for j = 1:length(lonall)
        knode = knode + 1;
        disp(['Processing node ', num2str(knode), ' of ', num2str(xnpts*ynpts)]);
        % Get velocity model at current location
        lat = latall(i);
        lon = lonall(j);
        % note that                 m0 starts at 0 km depth. The elevation infomration is
        % missing, to include topography set if_topo to 1
        if_topo = 1;
        [m0,~] = obtain_crust1_v2(lat,lon,[],if_topo);
        % m0 = obtain_ustc_litho(lat,lon);
        % define the interface in the crust 1.0 model
        m_interface = [];
        for l = 1:size(m0,1)
            if l == 1
                m_interface(l,:) = m0(l,:);
            else
                m_interface(2*(l-1),:) = [m0(l,1) m0(l-1,2:4)];
                m_interface(2*(l-1)+1,:) = m0(l,:);
            end
        end
        dmoho(knode) = m0(end,1);

        % Combine with AK135 model below Moho
        dmax = 1000;
        keepz = z > dmoho(knode) & z < dmax;
        m_interface = [m_interface; z(keepz) vp(keepz) vs(keepz) rho(keepz)];
        EPS = 1e-6;
        ztemp = m_interface(:,1);
        idisc = find( ztemp(1:end-1) == ztemp(2:end) );
        ztemp(idisc) = ztemp(idisc) - EPS;
        zpos = 0:dz:zmax;
        vptemp = interp1( ztemp, m_interface(:,2), zpos(1:end-1)+0.5*dz, 'linear','extrap');
        vstemp = interp1( ztemp, m_interface(:,3), zpos(1:end-1)+0.5*dz, 'linear','extrap');
        Z = [Z; zpos(1:end-1)'];
        X = [X; ones(size(vptemp'))*lon];
        Y = [Y; ones(size(vptemp'))*lat];
        VP = [VP; vptemp'];
        VS = [VS; vstemp'];

        vpgrid(:,j,i) = vptemp;
        vsgrid(:,j,i) = vstemp;

    end
end

VP = smooth3(vpgrid,'gaussian',[5 5 5]);
VS = smooth3(vsgrid,'gaussian',[5 5 5]);


[X,Y,Z] = meshgrid(lonall,latall,zpos(1:end-1));

xslice = 109.1;   
yslice = 42;
zslice = [];
figure;
subplot(121)
slice(X,Y,Z,permute(VP,[3 2 1]),xslice,yslice,zslice)
set(gca,'ZDir','reverse');
shading interp

subplot(122)
slice(X,Y,Z,permute(VS,[3 2 1]),xslice,yslice,zslice)
set(gca,'ZDir','reverse');
shading interp
