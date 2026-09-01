function plotStations(DataStruct, demFile, external_station_file)
% PLOTSTATIONS - Using MATLAB Mapping Toolbox to plot station locations
%                on a background DEM (digital elevation model).
%
% Usage:
%   plotStations(DataStruct, demFile)
%
% Inputs:
%   DataStruct - struct array with .StationInfo.stla, .StationInfo.stlo, .StationInfo.sta
%   demFile    - path to a .mat file containing DEM data:
%                demLat, demLon, demZ
%
% Prerequisite: MATLAB Mapping Toolbox
%
% Example:
%   plotStations(DataStruct, 'myDEM.mat')

%% 1. Load DEM data
if nargin < 2 || isempty(demFile)
    demFile = 'dem_data.mat';  % 
    
end
if nargin < 3
    external_station_file = [];
end
S = load(demFile);
demLat = S.demLat;  % [M x 1]
demLon = S.demLon;  % [N x 1]
demZ   = S.demZ;    % [M x N], Z(i,j) 与 (demLat(i), demLon(j)) 对应

%% 2. 
lats = [];
lons = [];
staname = {};
for n = 1:length(DataStruct)
    lats(n)    = DataStruct(n).StationInfo.stla;
    lons(n)    = DataStruct(n).StationInfo.stlo;
    staname{n} = DataStruct(n).StationInfo.sta;
end

% 
[staname, idx] = unique(staname);
lats = lats(idx);
lons = lons(idx);


latMin = min(lats) - 0.25;  
latMax = max(lats) + 0.25;
lonMin = min(lons) - 0.25;
lonMax = max(lons) + 0.25;

%  latMin < latMax, lonMin < lonMax
if latMin == latMax, latMax = latMax + 0.5; end
if lonMin == lonMax, lonMax = lonMax + 0.5; end

latlim = [latMin, latMax];
lonlim = [lonMin, lonMax];

%% 4. (worldmap or axesm)
figure('Position',[10 10 1000 800],'Name','StationMap_MappingToolbox','Color','w');
worldmap(latlim, lonlim);
set(gca,'FontSize',10,'FontWeight','normal');

%% 5. dem
geoshow(demLat, demLon, demZ, ...
    'DisplayType','texturemap');

demcmap(demZ);   
colorbar; 
hold on;

%% 6. station

geoshow(lats, lons, ...
    'DisplayType','point', ...
    'Marker','^', 'MarkerFaceColor','r','MarkerEdgeColor','k', ...
    'MarkerSize',15);

for i = 1:length(lats)
    textm(lats(i), lons(i), staname{i}, ...
        'VerticalAlignment','bottom','HorizontalAlignment','right', ...
        'FontSize',12,'Color','k','FontWeight','normal');
end

%% 7.
title('Station Map with DEM', 'FontSize',12,'FontWeight','bold');
gridm('on');    framem('on');   mlabel on;  plabel on; 

%% 8. 
fault_files = dir('./visualization/faults/*txt');
for k = 1:length(fault_files)
    faults = read_faults_gmt(fullfile(fault_files(k).folder, fault_files(k).name));
    for l = 1:length(faults)
        fault = faults{l};
        % fault(:,1) -> lon, fault(:,2) -> lat
        geoshow(fault(:,2), fault(:,1),'DisplayType','line','LineWidth',2,'Color','k');
    end
end

if ~isempty(external_station_file)
    fid = fopen(external_station_file, 'r');

    networks = {};
    latitudes = [];
    longitudes = [];
    stations = {};

    while ~feof(fid)
        line = fgetl(fid); %
        if isempty(line) || line(1) == '#'
            continue;
        end
        data = strsplit(line, '|');
        networks = [networks; data{1}]; %
        latitudes = [latitudes; str2double(data{3})]; %
        longitudes = [longitudes; str2double(data{4})]; %
        stations = [stations; data{2}];
    end

    fclose(fid);

    unique_networks = unique(networks);
    num_networks = length(unique_networks);

    colors = lines(num_networks);
    symbol = '^';

    load coastlines;
    plotm(coastlat, coastlon, 'k');

    countries = shaperead('landareas.shp', 'UseGeoCoords', true);
    geoshow(countries, 'EdgeColor', 'k', 'FaceColor', 'none');

    for i = 1:num_networks
        network = unique_networks{i};
        idx = find(strcmp(networks, network));
        scatterm(latitudes(idx), longitudes(idx), 100, colors(i, :), symbol, 'filled','MarkerEdgeColor','k'); % 绘制站点
        for j = 1:length(stations(idx))
            textm(latitudes(idx(j)), longitudes(idx(j)), ...
                [network, ' ', stations{idx(j)}], ...
                'Color', colors(i, :), 'FontSize', 8, ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left');
        end
    end

    hold off;
end
end
