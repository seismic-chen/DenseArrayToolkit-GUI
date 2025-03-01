function varargout = EventMaper(varargin)
% EventMaper: a function used to draw locations of events on map
% It requires zmap functions to run. Input argument varargin must
% contain "handles" with retrieved events, event-parameters, and
% the index of selected events.
%
% USAGE: handle = EventMaper(handles) or handle = EventMaper(myEvents,
% parameters, index). Here, handles must have handles.myEvent, 
% handles.p.area (eg [-90 90 -180 180], handles.selection (eg [1 2 3 
% 4 5 7 9...])

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Aug, 2003

if nargin < 1
    msgbox('Need retrieved and selected events, event parameters...','Warning');
elseif nargin==1
    events=handles.myEvent;
    area=handles.p.area;
    selection=handles.selection;
elseif nargin >=3
    events=varargin{1};
    area=varargin{2};
    selection=varargin{3};
end
% get input arguments

if isempty(selection)
    beep;
    msgbox('Select events first...', 'Warning');
    return;
end

if ~isa(events, 'MatEvent')
    beep;
    msgbox('Input must be one object of MatEvent class','Warning');
    return;
end;

% load events which will be displayed in listbox
if isempty(events)
    beep;
    msgbox('Retrieve event first...','Warning');
    return;
end
% if have not selected events, please do it first
hf=figure;
m_proj('miller','lon',[area(3) area(4)],'lat',[area(1) area(2)]);
% set map range
m_coast('patch',[0.2 1 .8]);
m_grid('box','fancy','tickdir','in');
% draw map frame
set(gcf,'Name','Event Map','NumberTitle','off');

colormap(hsv);
c=colormap;
mags=[];        % initialize the string for legend
hold on;
for i=1:length(selection)
    inx=selection(i);       % get the index of selected event
    ein=events.preferredOrigin(inx);        % get the event information
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