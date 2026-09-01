function consistent_earthquakes = filter_earthquakes_by_azimuth(slon, slat, elon, elat, max_angle_diff)
    % FILTER_EARTHQUAKES_BY_AZIMUTH 
    %
    % input:
    %   slon - 
    %   slat - 
    %   elon - 
    %   elat - 
    %   max_angle_diff -
    %
    % output:
    %   consistent_earthquakes -

    station_data = [slon, slat];

    mean_station = mean(station_data);
    centered_station = station_data - mean_station;

    % PAC
    [coeff, ~, ~] = pca(centered_station);
    line_direction = coeff(:, 1); 

    % 
    line_azimuth = atan2d(line_direction(1), line_direction(2)); % degree
    if line_azimuth < 0
        line_azimuth = line_azimuth + 360; %  [0, 360] 
    end

    % 
    line_azimuth_opposite = mod(line_azimuth + 180, 360); 

    % ）
    earthquake_data = [elon, elat];
    num_earthquakes = size(earthquake_data, 1);
    earthquake_azimuth = zeros(num_earthquakes, 1);

    for i = 1:num_earthquakes
        earthquake_azimuth(i) = azimuth(mean_station(2), mean_station(1), elat(i), elon(i));
    end

    % 
    angle_diff1 = abs(earthquake_azimuth - line_azimuth);
    angle_diff1 = min(angle_diff1, 360 - angle_diff1); 

    angle_diff2 = abs(earthquake_azimuth - line_azimuth_opposite);
    angle_diff2 = min(angle_diff2, 360 - angle_diff2); 

    % 
    angle_diff = min(angle_diff1, angle_diff2);

    % 
    consistent_earthquakes = find(angle_diff <= max_angle_diff);
    % pring
    fprintf('Number of events: %d \n', length(consistent_earthquakes));

    % plot
%     figure;
%     plot(station_data(:, 1), station_data(:, 2), 'bo', 'DisplayName', '台站位置');
%     hold on;
%     plot(earthquake_data(:, 1), earthquake_data(:, 2), 'k.', 'DisplayName', '地震位置');
%     plot(earthquake_data(consistent_earthquakes, 1), earthquake_data(consistent_earthquakes, 2), 'rx', 'DisplayName', '与测线一致的地震');
%     xlabel('Lat (deg)');
%     ylabel('Lon (deg)');
%     title('Distribution of events');
%     legend show;
%     axis equal;
%     grid on;
%     hold off;

%   Options:
%   - 'lambertstd': Lambert Azimuthal Equal-Area
%   - 'eqdazim': Azimuthal Equidistant
%   Here we choose Lambert to preserve area, or eqdazim to preserve distance.

figure('Position',[100 100 700 700], 'Name','Global station–event map', 'Color','w');
lat0 = mean(station_data(:,2));
lon0 = mean(station_data(:,1)); 
ax = axesm('ortho', 'Origin',[lat0 lon0 0], ...
           'Frame','on', 'Grid','on', 'MeridianLabel','on', 'ParallelLabel','on');
axis off;tightmap;

load topo
topo_shift = topo(:, [181:360 1:180]);
nlat = size(topo_shift,1);
nlon = size(topo_shift,2);
lat  = linspace(-90, 90,  nlat);
lon  = linspace(-180,180, nlon);
[lonGrid, latGrid] = meshgrid(lon, lat);
hTopo = surfm(latGrid, lonGrid, topo_shift);
set(hTopo, 'EdgeColor','none', 'FaceAlpha',0.3);         
demcmap(turbo); clim([-8000 8000])

coast = load('coastlines.mat');
geoshow(coast.coastlat, coast.coastlon, 'DisplayType','line', 'Color','k', 'LineWidth',1);

hold on;
plotm(station_data(:,2), station_data(:,1), ...
      'b^', 'MarkerFaceColor','b', 'MarkerSize',8, 'DisplayName','Stations');

scatterm(earthquake_data(:,2), earthquake_data(:,1), ...
         500, 'pk', 'filled', 'DisplayName','Earthquakes');

scatterm(earthquake_data(consistent_earthquakes,2), ...
         earthquake_data(consistent_earthquakes,1), ...
         800, 'pr', 'filled', 'DisplayName','Profile events');


end


