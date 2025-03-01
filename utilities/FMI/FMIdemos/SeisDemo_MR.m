% FMI Seismogram demo for multiple retrieving
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

% enable iterator and retrieve seismograms in multiple times
mySeis = enable_iterator(mySeis);
mySeis = init_iterator(mySeis,10);
nIter = 0;
nSeis = 0;
while(1)
    nIter = nIter + 1;
    [mySeis, nRetrieve] = seis_retrieve(mySeis,'selected');   % retrieve seismograms from the DMC
    mySeis  = channel_info(mySeis);         % retrieve the station information from the DMC
    nSeis = nSeis + nRetrieve;
    % save seismograms into files
    ot = string2time(eventOrigin.originTime);
    ot(6) = round(ot(6));
    tempFileName = sprintf('%04d-%02d-%02dT%02d-%02d-%02d_%03d.mat', ot, nIter);
    seismogram = get(mySeis,'seismogram');
    seisAttr = get(mySeis,'seismogramAttr');
    seisChanInfo = get(mySeis,'seisChanInfo');   % get channel info from the obect and put them into a structure
    eval(['save ' tempFileName ' seismogram seisAttr seisChanInfo']);
    
    iterator = get(mySeis,'iterator');
    if iterator.nRemained == 0
        break;
    end;
end;

if nSeis>0;
    %  -------  plot a record section of the seismograms -------
    hold off
    scal=0.5;                 % maximum amplitude of each seismogram scaled to 'scal' degrees of epicentral distance
    delta=zeros(nSeis,1);     % epicentral distances
    azim=zeros(nSeis,1);      % azimuths
    backazim=zeros(nSeis,1);  % back azimuths
    max_amp=delta;            % maximum amplitude of each seismogram
    kk = 0;                   % index for all seismograms
    for ii = 1:nIter;         % loop over each iteration
        ot = string2time(eventOrigin.originTime);
        ot(6) = round(ot(6));
        tempFileName = sprintf('%04d-%02d-%02dT%02d-%02d-%02d_%03d.mat', ot, ii);
        load(tempFileName);
        for jj = 1:length(seismogram)   % loop over each seismogram in an iteration
            kk = kk+1;
            [delta(kk),azim(kk),backazim(kk)]=delaz(...
                eventOrigin.location.latitude,eventOrigin.location.longitude, ...
                seisChanInfo(jj).location.latitude, seisChanInfo(jj).location.longitude, 0); 
            if length(seismogram{jj})>0; 
                max_amp(kk) = (max(seismogram{jj}) - min(seismogram{jj}) ) /2; 
                tim = timedif(eventOrigin.originTime,seisAttr(jj).beginTime) + ...
                    [0:seisAttr(jj).nPoint-1]*seisAttr(jj).sampleIntv*.001;   % time vector
                plot(tim, delta(kk) + (seismogram{jj}-mean(seismogram{jj}))/max_amp(kk)*scal)
                text(tim(1),delta(kk),sprintf('%s', seisChanInfo(jj).stationCode), 'horizontalalign','right')
                hold on
            end;
        end;
        
        delete(tempFileName);
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
