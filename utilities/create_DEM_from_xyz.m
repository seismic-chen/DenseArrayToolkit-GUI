function create_DEM_from_xyz(filename)
if nargin < 1
    filename = '../visualization/earth_relief_15s.xyz';
%     filename = '/Users/yunfeng/30_40/research/Baiyan_Ebo/figures/gmtplots/station_map/earth_relief_30s.xyz';
end
dem_model = load(filename);

lon = dem_model(:,1);
lat = dem_model(:,2);
elev = dem_model(:,3);

nlon = length(unique(lon));
nlat = length(unique(lat));

demLon = reshape(lon,nlon,nlat)';
demLat = reshape(lat,nlon,nlat)';
demZ = reshape(elev,nlon,nlat)';

figure; imagesc(demLon(1,:),demLat(:,1),demZ); axis xy;

demLat = demLat;  % [M x 1]
demLon = demLon;  % [N x 1]
demZ = demZ;    % [M x N], Z(i,j) 与 (demLat(i), demLon(j)) 对应

%save Qaidam_DEM.mat demLat demLon demZ;
save Baiyanebo_DEM.mat demLat demLon demZ;
