
function extract_dem(sta_lons,sta_lats)

fdir = './visualization/extract_dem/';

% sta_lons = [-116,-110];
% sta_lats = [43,45];

marginy = (max(sta_lats) - min(sta_lats))*0.1; 
marginx = (max(sta_lons) - min(sta_lons))*0.1;

latlim = [min(sta_lats) - marginy, max(sta_lats) + marginy];
lonlim = [min(sta_lons) - marginx, max(sta_lons) + marginx];

X = readhgt([latlim(1),latlim(2),lonlim(1),lonlim(2)], 'srtm3','interp','outdir',fdir); 
% X = readhgt([latlim(1),latlim(2),lonlim(1),lonlim(2)], 'srtm3','interp','plot'); 

yy = X.lat(:);
xx = X.lon(:);

[demLon,demLat] = meshgrid(xx,yy);
demZ = X.z;

delete([fdir,'*.hgt']);

save('./visualization/dem_data.mat', 'demLat', 'demLon', 'demZ');