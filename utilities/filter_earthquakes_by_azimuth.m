function consistent_earthquakes = filter_earthquakes_by_azimuth(slon, slat, elon, elat, max_angle_diff)
    % FILTER_EARTHQUAKES_BY_AZIMUTH 筛选与测线方向一致的地震
    %
    % 输入参数:
    %   slon - 台站的经度 (列向量)
    %   slat - 台站的纬度 (列向量)
    %   elon - 地震的经度 (列向量)
    %   elat - 地震的纬度 (列向量)
    %   max_angle_diff - 允许的最大方位角差值 (单位: 度)
    %
    % 输出参数:
    %   consistent_earthquakes - 与测线方向一致的地震索引

    % 将台站位置合并为一个矩阵
    station_data = [slon, slat];

    % 对台站数据进行中心化
    mean_station = mean(station_data);
    centered_station = station_data - mean_station;

    % 使用PCA找到最佳拟合测线的方向
    [coeff, ~, ~] = pca(centered_station);
    line_direction = coeff(:, 1); % 第一个主成分方向

    % 计算测线的方位角（从正北方向顺时针计算）
    line_azimuth = atan2d(line_direction(1), line_direction(2)); % 单位: 度
    if line_azimuth < 0
        line_azimuth = line_azimuth + 360; % 将方位角转换到 [0, 360] 范围
    end

    % 测线有两个方向，相差180度
    line_azimuth_opposite = mod(line_azimuth + 180, 360); % 反方向

    % 计算每个地震相对于测线中心点的方位角（考虑球面几何）
    earthquake_data = [elon, elat];
    num_earthquakes = size(earthquake_data, 1);
    earthquake_azimuth = zeros(num_earthquakes, 1);

    for i = 1:num_earthquakes
        % 使用球面几何计算方位角
        earthquake_azimuth(i) = azimuth(mean_station(2), mean_station(1), elat(i), elon(i));
    end

    % 计算地震方位角与测线方位角的差值（考虑两个方向）
    angle_diff1 = abs(earthquake_azimuth - line_azimuth);
    angle_diff1 = min(angle_diff1, 360 - angle_diff1); % 考虑圆周对称性

    angle_diff2 = abs(earthquake_azimuth - line_azimuth_opposite);
    angle_diff2 = min(angle_diff2, 360 - angle_diff2); % 考虑圆周对称性

    % 取两个方向中的最小差值
    angle_diff = min(angle_diff1, angle_diff2);

    % 筛选出与测线方向一致的地震
    consistent_earthquakes = find(angle_diff <= max_angle_diff);

    % 输出筛选结果
    fprintf('找到 %d 个与测线方向一致的地震。\n', length(consistent_earthquakes));

    % 可视化结果
%     figure;
%     plot(station_data(:, 1), station_data(:, 2), 'bo', 'DisplayName', '台站位置');
%     hold on;
%     plot(earthquake_data(:, 1), earthquake_data(:, 2), 'k.', 'DisplayName', '地震位置');
%     plot(earthquake_data(consistent_earthquakes, 1), earthquake_data(consistent_earthquakes, 2), 'rx', 'DisplayName', '与测线一致的地震');
%     xlabel('经度');
%     ylabel('纬度');
%     title('地震位置及其与测线方向的一致性');
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
set(hTopo, 'EdgeColor','none', 'FaceAlpha',0.3);         % 背景透明度 (0~1)
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


