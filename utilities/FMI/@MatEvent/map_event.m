function map_event(myEvent,mapTitle)

% MAP_EVENT plot events on map
% 
% In the map, events are depicted by filled circles. The size of circles
% represents magnitude. The color of the circles represents depth.
%
% Maps are drawn using m_map, matlab code written and freely distributed 
% by Richard Pawlowicz (http://www.eos.ubc.ca:/~rich).
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

if myEvent.nEvent==0
    beep;
    disp('  No event found!');
    return;
end;
if length(myEvent.name) == 0
    beep;
    disp('  No event information found! ');
end;

if ~exist('mapTitle')
    mapTitle=[];
end;

min_lat=9999;
max_lat=-9999;
min_lon=9999;
max_lon=-9999;
for ii=1:myEvent.nEvent
    lat=myEvent.preferredOrigin(ii).location.latitude;
    lon=myEvent.preferredOrigin(ii).location.longitude;
    if lat<min_lat min_lat=lat;end;
    if lat>max_lat max_lat=lat;end;
    if lon<min_lon min_lon=lon;end;
    if lon>max_lon max_lon=lon;end;
end;

diff_lon=max_lon-min_lon;
if diff_lon==0 diff_lon=1;end;
if diff_lon>180
    lonLimit=[-180 180];
else
    lonLimit=[min_lon-diff_lon*0.2 max_lon+diff_lon*0.2];
end;
diff_lat=max_lat-min_lat;
if diff_lat==0 diff_lat=1;end;
if diff_lon>180
    latLimit=[-90 90];
else
    latLimit=[min_lat-diff_lat*0.2 max_lat+diff_lat*0.2];
end;
selection=[1:myEvent.nEvent];
area=[latLimit lonLimit];

selection = 1:myEvent.nEvent;

% if have not selected events, please do it first
hf=figure(gcf);  %open current figure
clf
m_proj('miller','lon',[area(3) area(4)],'lat',[area(1) area(2)]);
% set map range
%m_coast('patch',[0.2 1 .8]);
if abs(diff(lonLimit))>90
    m_coast('patch',[0.2 1 .8]);
    %m_gshhs_c('patch',[0.2 1 .8]);
elseif abs(diff(lonLimit))>30
    m_gshhs_l('patch',[0.2 1 .8]);
elseif abs(diff(lonLimit))>5
    m_gshhs_i('patch',[0.2 1 .8]);
else
    m_gshhs_h('patch',[0.2 1 .8]);
end;

m_grid('box','fancy','tickdir','in');
if ~isempty(mapTitle)
    title(mapTitle);
end;
% draw map frame
set(gcf,'Name','Event Map','NumberTitle','off');

colormap(hsv);
c=colormap;
mags=[];        % initialize the string for legend
hold on;
for i=1:length(selection)
    inx=selection(i);       % get the index of selected event
    ein=myEvent.preferredOrigin(inx);        % get the event information
    lon = ein.location.longitude;
    lat = ein.location.latitude;
    cinx=c(max(1,round(2.5*sqrt(ein.location.depth))),:);
    mag = ein.magnitude.value(1);
    if mag<2 mag=2;end;
    if mag>8 mag=8;end;
 
    mags{i} = num2str(mag);
    H=m_plot(lon,lat,'o','markersize',round(mag*2),'markeredgecolor',cinx, ...
        'markerfacecolor',cinx);
    T=m_text(lon,lat,ein.originTime(1:10),'color',[0.5 0.1 0.2],'FontSize',8,'FontWeight','bold','VerticalAlignment','bottom');
end

a1=axes('position',[0.02 0.14 0.02 0.75]);
image([1:64]');
set(gca,'YAxisLocation','right');
set(gca,'YTick',[1 9 17 26 36 48 62]);
set(gca,'YTickLabel',{'0','10','40','100','200','350','600'});
set(gca,'XTick',[]);
title('Depth');
xlabel('km');

a2=axes('position',[0 0 1 0.14]);
plot(0.3,0.2,'o','markersize',4);text(0.31,0.2,'2.0-');
hold on
plot(0.4,0.2,'o','markersize',6);text(0.41,0.2,'3.0');
plot(0.5,0.2,'o','markersize',8);text(0.51,0.2,'4.0');
plot(0.6,0.2,'o','markersize',10);text(0.615,0.2,'5.0');
plot(0.7,0.2,'o','markersize',12);text(0.72,0.2,'6.0');
plot(0.8,0.2,'o','markersize',14);text(0.825,0.2,'7.0');
plot(0.9,0.2,'o','markersize',16);text(0.93,0.2,'8.0+ MB');
xlim([0 1.2]);
ylim([0 1]);
set(gca,'visible','off')
