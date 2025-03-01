function mySeis = channel_location(mySeis)

% CHANNEL_INFO retrieve channel locations only for avialable seismograms.
% Method seis_query must be called before call channel_info
%
% mySeis = channel_info(mySeis)
%

import edu.iris.Fissures.IfNetwork.*;
import edu.iris.Fissures.network.*;
import edu.iris.Fissures.*;

if isempty(mySeis.availSeisFilt)
    warning('Request filters for available seismogram is empty. Has seis_query been called first?');
end;

disp('Retrieving station locations ... ');

% sort availSeisFilt by station, site and channel name
rf = mySeis.availSeisFilt;
chanName = {};
for ii=1:length(rf)
    name=[];
    name=[char(rf(ii).channel_id.network_id.network_code) '.'];
    name=[name char(rf(ii).channel_id.station_code) '.'];
    name=[name char(rf(ii).channel_id.site_code) '.'];
    name=[name char(rf(ii).channel_id.channel_code)];
    chanName{ii} = name;   % station name
%     name=[name '.' mySeis.availFrom{ii}];
%     chanNameAll{ii} = name;  % station name with POND or SPYDER tag.
end;
% [chanNameAll, index] = sort(chanNameAll);
% chanName = chanName(index);

[chanName, index] = sort(chanName);
mySeis.availSeisFilt = mySeis.availSeisFilt(index);
% mySeis.availFrom = mySeis.availFrom(index);
rf = mySeis.availSeisFilt;

lastStation = [];
badChannel = [];
lastBad = 0;
for ii=1:length(rf)
    seisBeginTime = rf(ii).start_time;
    chanId = rf(ii).channel_id;
    chanId.begin_time = seisBeginTime;
    chanId.network_id.begin_time = seisBeginTime;
    
    if strcmp(chanId.station_code,lastStation) && (ii-1)~=lastBad
        availChanLoc(ii) = availChanLoc(ii-1);
    else
        try
            chan(ii) = mySeis.netExp.retrieve_channel(chanId);
            availChanLoc(ii).name = chanName{ii};
            availChanLoc(ii).longitude = chan(ii).my_site.my_station.my_location.longitude;
            availChanLoc(ii).latitude = chan(ii).my_site.my_station.my_location.latitude;
            availChanLoc(ii).elevation = chan(ii).my_site.my_station.my_location.elevation.value;
            availChanLoc(ii).elevationUnit = char(chan(ii).my_site.my_station.my_location.elevation.the_units.toString);
        catch
            disp(lasterr);
            badChannel = [badChannel ii];
            lastBad = ii;
            availChanLoc(ii).name = [];
            availChanLoc(ii).longitude = [];
            availChanLoc(ii).latitude = [];
            availChanLoc(ii).elevation = [];
            availChanLoc(ii).elevationUnit = [];
        end;   
    end
    lastStation = chanId.station_code;
    fprintf('.');
    if mod(ii,50) == 0
        fprintf(' %%%d\n',round(ii/length(rf)*100));
    end;
end;
disp('  Done!');

if ~isempty(badChannel)
    disp(['Warning: ' num2str(length(badChannel)) ' bad seismograms have been removed!']);
    availChanLoc(badChannel) = [];
    index = 1:length(mySeis.availSeisFilt);
    index(badChannel) = [];
    mySeis.availSeisFilt = mySeis.availSeisFilt(index);
end;

mySeis.availChanLoc = availChanLoc;
mySeis.selectedSeis = [];
