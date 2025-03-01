function varargout = SeisMaper(varargin)
% SeisMaper: a function used to draw locations of seismograms on map
% It requires zmap functions to run. Input arguments are paired with
% property/value.
% 
% Properties/Values
%       'Event', event. event is structural array containing:
%
%           time: a string for event originTime, e.g. '2001-01-09T16:49:28.000Z'
%           longitude: a numerical number for event longitude, e.g. 167.1700
%           latitude: a numerical number for event latitude, e.g. -14.9280
%           depth: a numerical number for event depth, e.g. 103
%
%       'Seisatt', seisatt. seisatt contains attributes of retrieved seismograms
%
%       'Channelinfo', channelinfo. channelinfo is about channel
%       information
%
%       'Selection', selection. selection is array of index of selected seismograms
%        (eg [1 2 3 4 5 7 9...])

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Aug, 2003
event=[];
seis=[];
selection=[];
% initialize

for i=1:nargin
    ts=varargin{i};
    if ischar(ts)
        switch lower(ts)
            case 'event'
                event=varargin{i+1};
            case 'seisatt'
                seisatt=varargin{i+1};
            case 'channelinfo'
                channelinfo=varargin{i+1};
            case 'selection'
                selection=varargin{i+1};
        end;
    end;
end;
% get input arguments
 
if isempty(selection)
    selection=[1:length(seisatt)];
end

area=[-90 90 -180 180];
chan=channelinfo;
att=seisatt;
% if have not selected events, please do it first
hf=figure;
m_proj('miller','lon',[area(3) area(4)],'lat',[area(1) area(2)]);
% set map range
m_coast('patch',[0.2 1 .8]);
m_grid('box','fancy','tickdir','in');
% draw map frame
set(gcf,'Name','Seismogram Map','NumberTitle','off');

colormap(hsv);
c=colormap;
hold on;

if exist('event') && ~isempty(event)
    cinx=c(max(1,round(2.5*sqrt(event.depth))),:);
    H=m_plot(event.longitude,event.latitude,'o','markersize',16,'markeredgecolor',cinx, ...
        'markerfacecolor',cinx);
    T=m_text(event.longitude,event.latitude,{event.time(1:10);event.time(12:23)},...
        'color',[0 0.2 0.9],'FontSize',8,'FontWeight','bold','VerticalAlignment','bottom');
end

blon=-400;
blat=-400;
% backup for longitude and latitude in order check duplicates

for i=1:length(selection)
    inx=selection(i);       % get the index of selected seismograms
    lon = chan(inx).location.longitude;
    lat = chan(inx).location.latitude;
    if lat==blat && lon==blon
        continue
    end
    
    ss=[];
    if ~isempty(event)
        stx=timedif(event.time,att(inx).beginTime);
        ss=['T:',num2str(round(stx))];
        
        [delta,azim,backazim]=delaz(event.latitude,event.longitude, ...
            chan(inx).location.latitude,  chan(inx).location.longitude, 0);
        ss=[ss,' D:',num2str(round(delta)),' A:',num2str(round(azim))];
    end
    
    ss={ss; ['ID-',att(inx).network,'.',att(inx).station,'.',att(inx).site]}; %,'.',att(inx).channel]};
    
    H=m_plot(lon,lat,'o','markersize',6,'markeredgecolor',[0.9 0.1 0.2], ...
        'markerfacecolor',[0.9 0.1 0.2]);
    T=m_text(lon,lat,ss,'color',[0.9 0.1 0.2],'FontSize',8,'FontWeight','bold','VerticalAlignment','bottom');
    blon=lon;
    blat=lat;
end

if ~isempty(event)
    a1=axes('position',[0.02 0.14 0.02 0.75]);
    image([1:64]');
    set(gca,'YAxisLocation','right');
    set(gca,'YTick',[1 9 17 26 36 48 62]);
    set(gca,'YTickLabel',{'0','10','40','100','200','350','600'});
    set(gca,'XTick',[]);
    title('Depth');
    xlabel('km');
end