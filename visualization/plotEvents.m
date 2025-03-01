function plotEvents(DataStruct)
% PLOTEVENTS - Plot seismic events in the 30-90 degree range from a chosen station,
%              using MATLAB Mapping Toolbox with an azimuthal projection.
%
% Usage:
%   plotEvents(DataStruct)
%
% DataStruct needs fields:
%   DataStruct(i).EventInfo.evla, evlo, evid    % event lat, lon, ID
%   [Optionally a chosen station lat/lon, see below]
%
% This example:
%   1) Takes a 'reference station' lat/lon (hardcoded or from DataStruct),
%   2) Filters events whose great-circle distance is in [30°, 90°],
%   3) Uses an azimuthal projection (Lambert or Equidistant) with that station as origin,
%   4) Plots events in the radial coordinate.

% -------------------------------------------------------------------------
%% 1. Extract event lat/lon and ID from DataStruct
lats = [];
lons = [];
evid = {};
for n = 1:length(DataStruct)
    lats(n)    = DataStruct(n).EventInfo.evla;
    lons(n)    = DataStruct(n).EventInfo.evlo;
    if isfield(DataStruct(n).EventInfo,'evid')
        evid{n} = DataStruct(n).EventInfo.evid;
    else
        evid{n} = num2str(n);  % fallback if evid not provided
    end
end

% 去重（若有相同ID的事件）
[evid, idx] = unique(evid);
lats = lats(idx);
lons = lons(idx);

% -------------------------------------------------------------------------
%% 2. Define the reference station
% 你需要根据实际情况获取台站经纬度:
% 1) 如果 DataStruct 中有 StationInfo 的“中心台站”，可从那里获取
% 2) 如果只想示例，就手动指定一个stationLat, stationLon
%   (例如 stationLat=34, stationLon=103)

% [示例] 假设 DataStruct(1) 里有一个中心台信息:
stationLat = DataStruct(1).StationInfo.stla;
stationLon = DataStruct(1).StationInfo.stlo;

% -------------------------------------------------------------------------
%% 3. Compute the great-circle distance from the station to each event
%   distance(lat1, lon1, lat2, lon2, spheroid)
%   - lat/lon in degrees, returns distance in degrees by default
%   - You can specify 'degrees' or pass a referenceEllipsoid if needed
distDeg = distance(stationLat, stationLon, lats, lons);

% -------------------------------------------------------------------------
%% 4. Filter out events outside [30°, 90°]
idxRange = (distDeg >= 30 & distDeg <= 90);

if ~any(idxRange)
    warning('No events found in [30°, 90°] range from the station.');
    % 也可选择直接return
end

latsSel = lats(idxRange);
lonsSel = lons(idxRange);
distSel = distDeg(idxRange);
evidSel = evid(idxRange);

% -------------------------------------------------------------------------
%% 5. Create a figure & set up an azimuthal map projection
%   Options:
%   - 'lambertstd': Lambert Azimuthal Equal-Area
%   - 'eqdazim': Azimuthal Equidistant
%   Here we choose Lambert to preserve area, or eqdazim to preserve distance.

figure('Position',[10 10 800 800],'Name','Teleseismic Events','Color','w');

% Create map axes with Lambert Azimuthal (good for polar or local region).
% 'origin' sets the center of projection as [Lat, Lon]
%   The first two elements of origin are the lat & lon of the center,
%   the third element can be the rotation, usually 0.
% axesm('MapProjection','lambertstd','MapLatLimit',[0 90],...
%       'Origin',[stationLat stationLon 0]);


% 若想用等距方位投影:
axesm('eqdazim', 'Origin',[stationLat stationLon 0]);

% MapLatLimit = [0 90] 表示半径90度范围(从中心到 90°)。 
% 你也可以写 [30 90] 如果只想显示内环-外环之间, 
% 但会导致图上不显示 <30° 区域. 看需求决定.

% 检查当前坐标轴是否是地图坐标轴
if isempty(gcm)
    error('Current axes is not a map axes.');
end

hold on;

% 打开格网、边框、标签
gridm('on');   % 绘制经纬网
mlabel('on');  % 纬度标注
plabel('on');  % 经度标注
framem('on');  % 地图框

title(sprintf('Events within [30°, 90°] Range'), ...
      'FontSize',12,'FontWeight','bold');
% -------------------------------------------------------------------------
%% 6. Plot the center station as a special marker
%   You can use plotm(stationLat, stationLon, markerSpec) or geoshow.
plotm(stationLat, stationLon, '^','MarkerSize',12,'MarkerFaceColor','b');
textm(stationLat, stationLon, '  Station','FontSize',10,'Color','b');
pause(0.05)
% -------------------------------------------------------------------------
%% 7. Plot the selected events
%   We can use scatterm or plotm
scatterm(latsSel, lonsSel, 200, 'r', 'p','filled', 'MarkerEdgeColor','k');
% 这里 color 由 distSel 决定, 你也可以固定颜色

% 标注事件 ID (可选)
for i = 1:length(latsSel)
    textm(latsSel(i), lonsSel(i), ['  ' evidSel{i}], ...
          'VerticalAlignment','middle','FontSize',8,'Color','k');
end

% 颜色栏(可选) - 显示距离
colormap(jet);
caxis([30, 90]);
% colorbar('Ticks',[30,50,70,90],'TickLabels',{'30°','50°','70°','90°'}, ...
%          'FontSize',8,'Location','eastoutside',...
%          'LabelString','Distance (deg)');

%% ============ 8. 绘制海岸线 =============
% 1) 可以使用内置 coastlines.mat
coast = load('coastlines.mat');  % coast.coastlat, coast.coastlon
% 或 2) geoshow('landareas.shp','FaceColor',[0.7 0.9 0.5]);
geoshow(coast.coastlat, coast.coastlon, ...
        'DisplayType','line','Color','k','LineWidth',1);
% -------------------------------------------------------------------------
% Done!
end