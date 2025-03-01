function profileStruct=getImagingProfile(DataStruct,profile_length)
if nargin < 2 
    profile_length = 4;
end
profileStruct = struct('line_points', [], ...
    'projected_points', [], ...
    'center', [], ...
    'direction', [], ...
    'rx', []) ;

% 定义观测系统
stationList = getStations(DataStruct);
stlo = [stationList.stlo]';
stla = [stationList.stla]';

% project the stations onto a profile
[line_points, projected_points, center, direction] = fit_and_project_stations(stlo,stla,profile_length,1);
slatp=projected_points(:,2);
slonp=projected_points(:,1);
lat1 = line_points(1,2);
lon1 = line_points(1,1);

[deg0,~]= distance(lat1,lon1,slatp,slonp);
rx = deg0*2*pi*6371/360;

profileStruct.line_points = line_points;
profileStruct.projected_points = projected_points;
profileStruct.center = center;
profileStruct.direction = direction;
profileStruct.rx = rx;