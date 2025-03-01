function [mySeis, nSelect] = seis_select(mySeis, origin, chanArea, timeWindow , twin,  handles)

% SEIS_SELECT selected seismograms by channel location, time window and
% database.
%
% [mySeis nSelect] = seis_select(mySeis, origin, chanArea, timeWindow, catalog)
% 
% Input arguments
%    origin:  Earthquake origin. This is a structure defined as follow
%             origin.time       (string or FMI time structure)
%                   .longitude  (in degree)
%                   .latitude   (in degree)
%                   .depth      (in km)
%
%  chanArea:  The interested area of channel location. The area can be
%             either specified by a square box or distance-azimuth. This is
%             a structure defined as follow
%             chanArea.type     ('boxArea' or 'distAzim')
%                     .minLon   (if type is 'boxArea')
%                     .maxLon   (if type is 'boxArea')
%                     .minLat   (if type is 'boxArea')
%                     .maxLat   (if type is 'boxArea')
%                     .minDist  (if type is 'distAzim')
%                     .maxDist  (if type is 'distAzim')
%                     .minAzim  (if type is 'distAzim')
%                     .maxAzim  (if type is 'distAzim')
% 
%  timeWindow:  Only retrieve data for interesting time window. The time
%               window is specified by seismic phase and relative time to
%               arrival time of that phase. This is structure defined as follow
%               timeWindow.fromPhase   (string, single or multiple phase names)
%               timeWindow.fromRelTim  (in seconds)
%               timeWindow.toPhase     (string, single or multiple phase names)
%               timeWindow.toRelTime   (in seconds)
%
% This methods will only return total number of select seismogram request
% filters. The index of selected filters will save into class field
% 'selectedSeis'. Option 'chanArea' and 'timeWindow' may be empty for no
% selection on area or time window.
%
% Note that time in arguments could be either string or FMI time structure.
%
% Phase Naming: Since this method calls MatTauP/TauP toolkit to calculate 
% travel time of seismic phases. All phase names must compliant with naming
% rules of Taup toolkit. Use 'help phaseName' to look for details. Please  
% note that PKPab, PKPcd, PKPdf are not acceptable. Instead, both PKPab and
% PKPdf are called as PKP; PKPdf is called PKIKP. So when PKP is specified,
% there might be multiple arrivals. To cover all distance range, we suggest 
% to use phase sequence 'p,P,Pdiff,PKIKP' for first P arrival, and use phase
% sequence 's,S,SKS,SKIKS' for first S-wave arrival. If there are multiple 
% arrivals, the earliest arrival will be used for calculating start of the 
% time window and the later arrival will be used for calculating end of the 
% time window. 
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also phaseName
%


% check input arguments
if nargin<2
    error('At least three arguments are required');
end;
if exist('catalog')~=1 || isempty(catalog)
    catalog = 'farm_priority';
end;
if exist('chanArea')~=1 || isempty(chanArea)
    chanArea = [];
end;
if exist('timeWindow')~=1 || isempty(timeWindow)
    timeWindow = [];
end;

disp('Selecting seismograms ...');

areaSelect = ones(1, length(mySeis.availSeisFilt));
timeSelect = ones(1, length(mySeis.availSeisFilt));
%fromSeletc = ones(1, length(mySeis.availSeisFilt));

if isempty(mySeis.availChanLoc)
    error('Please retrieve channel locations first');
end;
chanLon = [mySeis.availChanLoc.longitude];
chanLat = [mySeis.availChanLoc.latitude];

% select by channel locations
if ~isempty(chanArea)
    switch lower(chanArea.type)
        case 'boxarea'
            minLon = chanArea.minLon;
            if minLon>180 minLon=minLon-360;end;
            maxLon = chanArea.maxLon;
            if maxLon>180 maxLon=maxLon-360;end;
            if minLon>maxLon
                temp = minLon;
                minLon = maxLon;
                maxLon = temp;
            end;
            minLat = chanArea.minLat;
            maxLat = chanArea.maxLat;
            if minLat>maxLat
                temp = minLat;
                minLat = maxLat;
                maxLat = temp;
            end;
            
            areaSelect = (chanLon>=minLon & chanLon<=maxLon & chanLat>=minLat & chanLat<=maxLat);
        case 'distazim'
            minDist = chanArea.minDist;
            maxDist = chanArea.maxDist;
            if minDist>maxDist
                temp = minDist;
                minDist = maxDist;
                maxDist = temp;
            end;
            minAzim = chanArea.minAzim;
            maxAzim = chanArea.maxAzim;
            if minAzim>maxAzim
                temp = minAzim;
                minAzim = maxAzim;
                maxAzim = temp;
            end;
            
            [distance, azimuth, backAzim] = delaz(origin.latitude, origin.longitude, chanLat, chanLon, 0);
            areaSelect = (distance>=minDist & distance<=maxDist & azimuth>=minAzim & azimuth<=maxAzim);
        otherwise
            error('Wrong area select type');
    end;
end;
    
select = find(areaSelect);
for jj = 1:length(select)
        ii = select(jj);
        mySeis.availSeisFilt(ii).start_time = edu.iris.Fissures.Time(time2string(twin.fromRelTime),0);
        mySeis.availSeisFilt(ii).end_time = edu.iris.Fissures.Time(time2string(twin.toRelTime),0);
end

% check phase selection radiobutton is set
if get(handles.phase_select,'Value') && ~isempty(select) && ~isempty(timeWindow)
disp('  Calculate travel times for each channel ...');
    lastStaName = [];
    lastSelect = 0;
    for jj = 1:length(select)
        ii = select(jj);
        fprintf('.');
        if mod(jj,50) == 0
            fprintf(' %2d%%\n', round(jj/length(select)*100));
        end;

        dotIndex = findstr('.',mySeis.availChanLoc(ii).name);
        staName = mySeis.availChanLoc(ii).name(dotIndex(1)+1:dotIndex(2)-1);
        if strcmp(staName, lastStaName)
            timeSelect(ii) = lastSelect;
            %disp(['skipped  ', num2str(jj)]);
            continue;
        else
            lastStaName = staName;
            %disp(jj);
        end;
        
        fromTT = taupTime('iasp91', origin.depth, timeWindow.fromPhase, 'station', [chanLat(ii) chanLon(ii)], ...
            'event', [origin.latitude, origin.longitude]);
        if isempty(fromTT)
            warning('Phase does not exist at this distance or phase name is not correct !');
            fprintf('  Phase: %s,  channel: %s \n',timeWindow.fromPhase, mySeis.availChanLoc(ii).name);
            timeSelect(ii) = 0;
            lastSelect = timeSelect(ii);
            continue;
        end;
        ts = [];
        for k = 1:length(fromTT)
            ts(k) = fromTT(k).time;
        end;
        desireWinStart = timeshift(origin.time, min(ts)+timeWindow.fromRelTime);
        
        toTT = taupTime('iasp91', origin.depth, timeWindow.toPhase, 'station', [chanLat(ii) chanLon(ii)], ...
            'event', [origin.latitude, origin.longitude]);
        if isempty(toTT)
            %warning('Phase does not exist at this distance or phase name is not correct !');
            fprintf('  Phase: %s,  channel: %s \n',timeWindow.toPhase, mySeis.availChanLoc.name);
            timeSelect(ii) = 0;
            lastSelect = timeSelect(ii);
            continue;
        end;
        ts = [];
        for k = 1:length(toTT)
            ts(k) = toTT(k).time;
        end;
        desireWinEnd = timeshift(origin.time, max(ts)+timeWindow.toRelTime);
        
        availWinStart = char(mySeis.availSeisFilt(ii).start_time.date_time);
        availWinEnd = char(mySeis.availSeisFilt(ii).end_time.date_time);
        if timecmp(desireWinStart, availWinStart) > 0 % desireWinStart is later
            actualWinStart = desireWinStart;
        else
            actualWinStart = availWinStart;
        end;
        if timecmp(desireWinEnd, availWinEnd) < 0 % desireWinEnd is earlier
            actualWinEnd = desireWinEnd;
        else
            actualWinEnd = availWinEnd;
        end;
        
        if timedif(actualWinStart, actualWinEnd)<1 % the time window is short than 1 second
            timeSelect(ii) = 0;  % not selected
        else
            timeSelect(ii) = 1;
            mySeis.availSeisFilt(ii).start_time = edu.iris.Fissures.Time(time2string(actualWinStart),0);
            mySeis.availSeisFilt(ii).end_time = edu.iris.Fissures.Time(time2string(actualWinEnd),0);
        end;
        
        lastSelect = timeSelect(ii);
    end;
    disp('  Done!');
end;

% mySeis.selectedSeis = find(areaSelect & timeSelect & fromSelect);

mySeis.selectedSeis = find(areaSelect & timeSelect);
nSelect = length(mySeis.selectedSeis);

% reset retrieve iterator
mySeis = init_iterator(mySeis, 0);

fprintf('  %d seismograms selected. Done! \n',nSelect);
