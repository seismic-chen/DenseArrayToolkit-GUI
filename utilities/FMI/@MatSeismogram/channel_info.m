function mySeis = channel_info(mySeis)

% CHANNEL_INFO retrieve channel information of retrieved seismograms.
% Method seis_retrieve must be called before call channel_info
%
% mySeis = channel_info(mySeis) 
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%


import edu.iris.Fissures.IfNetwork.*;
import edu.iris.Fissures.network.*;
import edu.iris.Fissures.*;

if isempty(mySeis.seismogram)
    warning('Seismograms have not been retrieved');
    return;
end;

disp('Retrieving channel information ... ');

lastStation = [];
badChannel = [];
for ii=1:length(mySeis.seismogram)
    try
        chanId = mySeis.seisChanId(ii);
        chan(ii) = mySeis.netExp.retrieve_channel(chanId);
        
        chanInfo(ii).name = char(chan(ii).name);
        chanInfo(ii).networkCode = char(chan(ii).get_id.network_id.network_code);
        chanInfo(ii).stationCode = char(chan(ii).get_id.station_code);
        chanInfo(ii).citeCode = char(chan(ii).get_id.site_code);
        chanInfo(ii).channelCode = char(chan(ii).get_id.channel_code);
        chanInfo(ii).samplingIntv = chan(ii).sampling_info.interval.value/chan(ii).sampling_info.numPoints;
        chanInfo(ii).samplingIntvUnit = char(chan(ii).sampling_info.interval.the_units.toString);
        location.longitude = chan(ii).my_site.my_station.my_location.longitude;
        location.latitude = chan(ii).my_site.my_station.my_location.latitude;
        location.elevation = chan(ii).my_site.my_station.my_location.elevation.value;
        location.elevationUnit = char(chan(ii).my_site.my_station.my_location.elevation.the_units.toString);
        chanInfo(ii).location = location;
        orientation.azimuth = chan(ii).an_orientation.azimuth;
        orientation.dip = chan(ii).an_orientation.dip;
        chanInfo(ii).orientation = orientation;
        chanInfo(ii).startTime = char(chan(ii).effective_time.start_time.date_time);
        chanInf0(ii).endTiem = char(chan(ii).effective_time.end_time.date_time);
    catch
        disp(lasterr);
        badChannel = [badChannel ii];
        chanInfo(ii).name = [];
        chanInfo(ii).networkCode = [];
        chanInfo(ii).stationCode = [];
        chanInfo(ii).citeCode = [];
        chanInfo(ii).channelCode = [];
        chanInfo(ii).samplingIntv = [];
        chanInfo(ii).samplingIntvUnit = [];
        chanInfo(ii).location = [];
        chanInfo(ii).orientation = [];
        chanInfo(ii).startTime = [];
        chanInf0(ii).endTiem = [];
    end;
    
    fprintf('.');
    if mod(ii,50) == 0
        fprintf(' %%%d\n',round(ii/length(mySeis.seismogram)*100));
    end;
end;
disp('  Done!');

if ~isempty(badChannel)
    disp(['Warning: ' num2str(length(badChannel)) ' bad channels found!']);
end;

mySeis.seisChanInfo = chanInfo;

