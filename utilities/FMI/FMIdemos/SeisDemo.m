% FMI Seismogram demo
%
% Written by:
%   Qin Li and Ken Creager
%   Unverisity of Washington
%   qinli@u.washington.edu
%   kcc@ess.washington.edu
%   October, 2003
%

%---- Open Corba connection to DMC and initialize object of class MatEvent -------
clear myEvent mySeis
myEvent=MatEvent;

% ---- select big deep events from the POND for the year 2001
disp(' ');
disp('Searching for big deep events for the year 2001 in POND ...');
myEvent=find_event(myEvent,'Area',[-90 90 -180 180],'magnitude',[6 8], 'searchTypes',{'MB','MW'}, ... 
    'timeRange',{'2001-01-01T00:00:00.000Z','2001-12-28T23:59:59.000Z'},'depth',[300 700],'catalogs',{'POND'});
myEvent=retrieve_event_info(myEvent);    %  get the event information
list_events(myEvent);                    % list the event information to the screen

eventSelect = 4;
eventOrigin = myEvent.preferredOrigin(eventSelect);  % put the event information from event 4 into a structure

% check data availability (0-30 min relative to origin time)
timeWindow.startTime = string2time(eventOrigin.originTime);
timeWindow.endTime = timeshift(eventOrigin.originTime, 30*60);

% ------ initialize an object of class MatSeismogram ------
mySeis = MatSeismogram('Pond');

% - ask for seismograms from all avaialble networks starting with I, all stations, all site codes and BHZ 
mySeis = seis_query(mySeis,timeWindow,'I*','*','','BHZ');

%  --- get station latitude and longitude for all seismograms
mySeis = channel_location(mySeis);

origin.time      = eventOrigin.originTime;
origin.longitude = eventOrigin.location.longitude;
origin.latitude  = eventOrigin.location.latitude;
origin.depth     = eventOrigin.location.depth;
chanArea.type    = 'distAzim';
chanArea.minDist = 60;  %degree
chanArea.maxDist = 90;
chanArea.minAzim = 0;
chanArea.maxAzim = 360;
timeWindow.fromPhase   = 'p,P,Pdiff,PKIKP';  % first P
timeWindow.fromRelTime = -60;                % seconds
timeWindow.toPhase     = 's,S,SKS,SKIKS';    % first S 
timeWindow.toRelTime   = 1200;
%database               = 'Farm_Priority';    % take data from SPYDER only if not available from FARM
[mySeis nSelect] = seis_select(mySeis, origin, chanArea, timeWindow);

[mySeis, nRetrieve] = seis_retrieve(mySeis,'selected');   % retrieve seismograms from the DMC

if nRetrieve>0;
    mySeis       = channel_info(mySeis);         % retrieve the station information from the DMC
    seismogram   = get(mySeis,'seismogram');     % get seismograms from the object and put them into a cell array
    seisAttr     = get(mySeis,'seismogramAttr'); % get seismogram attribures from the object and put them into a structure
    seisChanInfo = get(mySeis,'seisChanInfo');   % get channel info from the obect and put them into a structure
    nseis=length(seismogram);
    
    %  -------  plot a record section of the seismograms -------
    hold off
    scal=0.5;                 % maximum amplitude of each seismogram scaled to 'scal' degrees of epicentral distance
    delta=zeros(nseis,1);     % epicentral distances
    azim=zeros(nseis,1);      % azimuths
    backazim=zeros(nseis,1);  % back azimuths
    max_amp=delta;            % maximum amplitude of each seismogram
    for k=1:nseis;            % loop over each seismogram
        [delta(k),azim(k),backazim(k)]=delaz(...
            eventOrigin.location.latitude,eventOrigin.location.longitude, ...
            seisChanInfo(k).location.latitude, seisChanInfo(k).location.longitude, 0); 
        if length(seismogram{k})>0; 
            max_amp(k) = (max(seismogram{k}) - min(seismogram{k}) ) /2; 
            tim = timedif(eventOrigin.originTime,seisAttr(k).beginTime) + ...
                [0:seisAttr(k).nPoint-1]*seisAttr(k).sampleIntv*.001;   % time vector
            plot(tim, delta(k) + (seismogram{k}-mean(seismogram{k}))/max_amp(k)*scal)
            text(tim(1),delta(k),sprintf('%s', seisChanInfo(k).stationCode), 'horizontalalign','right')
            hold on
        end
    end
    % ---- calculate and plot travel time curves ----
    tt=taupCurve('ak135',eventOrigin.location.depth,'P,pP,sP,PP,PcP,S,sS,ScS,SKS');
    for k=1:length(tt); plot(tt(k).time,tt(k).distance,':r');end
    ylim([min(delta)-scal,max(delta)+scal]);
    xlim([-Inf 1600])
    ylabel('Distance (deg)');
    xlabel('Time after origin time (s)')    
    event_title = sprintf('%s   %.3fN   %.3fE   %.1fkm',eventOrigin.originTime,...
        eventOrigin.location.latitude,eventOrigin.location.longitude,eventOrigin.location.depth);
    disp(event_title)
end
