function map_station(myNetwork,option)

% MAP_STATION plot events on map
%
% Usage
%   map_station(myNetwork,'networkStation') map stations in "networkStation"
%   map_station(myNetwork,'stationPool') map stations in "stationPool".
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

switch lower(option)
case 'networkstation'
    station=myNetwork.networkStation;
case 'stationpool'
    station=myNetwork.stationPool;
otherwise
    error('Unregnized option');
end;

if isempty(station)
    disp('No station found');
    return;
end;
min_lat=9999;
max_lat=-9999;
min_lon=9999;
max_lon=-9999;
for ii=1:length(station)
    lat=station(ii).location.latitude;
    if lat<min_lat min_lat=lat;end;
    if lat>max_lat max_lat=lat;end;
    lon=station(ii).location.longitude;
    if lon<min_lon min_lon=lon;end;
    if lon>max_lon max_lon=lon;end;
end;

diff_lon=max_lon-min_lon;
if diff_lon==0 diff_lon=1;end;
if diff_lon>180
    LonLimit=[-180 180];
else
    LonLimit=[min_lon-diff_lon*0.2 max_lon+diff_lon*0.2];
end;

diff_lat=max_lat-min_lat;
if diff_lat==0 diff_lat=1;end;
if diff_lon>180
    LatLimit=[-90 90];
else
    LatLimit=[min_lat-diff_lat*0.2 max_lat+diff_lat*0.2];
end;

xr=diff(LonLimit)*cos(pi/180*min(abs(floor(LatLimit(1)):ceil(LatLimit(2)))));
yr=diff(LatLimit);
ratio=xr/yr;
if ratio<1.5
    center=mean(LonLimit);
    range=diff(LonLimit);
    LonLimit(1)=center-0.5*range/ratio*1.5;
    LonLimit(2)=center+0.5*range/ratio*1.5;
    if LonLimit(1)<-180 LonLimit(1)=-180;end;
    if LonLimit(2)>180 LonLimit(2)=180;end;        
end;
if ratio>2
    center=mean(LatLimit);
    range=diff(LatLimit);
    LatLimit(1)=center-0.5*range*ratio/1.8;
    LatLimit(2)=center+0.5*range*ratio/1.8;
    if LatLimit(1)<-90 LatLimit(1)=-90;end;
    if LatLimit(2)>90 LatLimit(2)=90;end;
end;
xr=diff(LonLimit)*cos(pi/180*min(abs(floor(LatLimit(1)):ceil(LatLimit(2)))));
yr=diff(LatLimit);
    
h=figure(101);clf
h_height=450;
h_width=h_height*xr/yr;
if h_width>1000 h_width=1000;end;
set(h,'position',[20,700-h_height-20,h_width,h_height]);
if diff(LatLimit)<15 & diff(LonLimit)<15
    worldmap('hi',LatLimit,LonLimit);
else
    worldmap('lo',LatLimit,LonLimit);
end;
title('Station Map');

for ii=1:length(station)
    lat=station(ii).location.latitude;
    lon=station(ii).location.longitude;
    code=station(ii).code;
    plotm(lat,lon,'g^','markersize',6,'markerfacecolor',[0 1 0]);
    %textm(lat+diff(LatLimit)*0.02,lon+diff(LonLimit)*0.02,code,'color',[1 0 1],'FontSize',8);
    textm(lat,lon,code,'color',[1 0 1],'FontSize',8,'FontWeight','bold','VerticalAlignment','bottom');
end;
%hidem(gca);

return;


