function myNetwork=find_stationChannel(myNetwork, netId, channelOption,availableTime)

% FIND_STATIONCHANNEL find channels for stations in data member networkStation
%
% Usage: myNetwork=find_stationChannel(myNetwork,channel_option,availableTime)
%
% Input arguments:
%   "channelOption" must be string of length 3 to specify components. 
%   Asterisk could be used for substitution of any charactor
%   Example:
%       'BH*':  all components of broad band
%       'BHZ':  z component of broad band
%       '**Z':  z component of all period
%       'ALL':  all channels
%   "AvailableTime": If availableTime is specified, only the channels that
%   are effective at availableTime will be selected
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

import edu.iris.Fissures.IfNetwork.*;
import edu.iris.Fissures.network.*;
import edu.iris.Fissures.*;
import edu.iris.Fissures.model.*;
import edu.iris.Fissures.utility.*;
import java.io.*;
import java.lang.*;
import java.util.*;

if myNetwork.findStationDone~=1
    disp('Please run FIND_STATION first');
    return;
end;

% input arguments checking ...
if exist('availableTime')~=1
    availableTime = [];
end;

if nargin<2
    error('At least 2 input arguments needed');
end;
if ~isempty(availableTime) & ~isa(availableTime,'char') & ~isa(availableTime,'numeric')
    error('AvailableTime must be Fissures time string or time structure');
end;
if length(channelOption)~=3 | ~isa(channelOption,'char')
    error('channelOption must be 3 charactors string');
end;
channelOption=upper(channelOption);

% begin query
all_chanNum=0;
jNetId = convert_id(netId,'net');
netAccess = myNetwork.networkFnd.retrieve_by_id(jNetId);

for kk=1:myNetwork.networkStationNum
    fprintf('Retrieving channels for station %4s ... \n',myNetwork.networkStation(kk).id.stationCode);
    jStationId = convert_id(myNetwork.networkStation(kk).id,'sta');
    jChannels=netAccess.retrieve_for_station(jStationId);
    
%     if kk == 6
%         disp('debug');
%         for ii = 1:length(jChannels)
%             disp([int2str(ii) '   ' (jChannels(ii).name.toCharArray)']);
%         end;
%     end;

    cn=0;  % index of channel
    clear stationChannel
    for ii=1:length(jChannels)
        channel_chosen=0;
        myChannel.id = convert_id(jChannels(ii).get_id, 'chan');
        myChannel.fullCode = [(jChannels(ii).get_id.network_id.network_code.toCharArray)' '.' ...
                (jChannels(ii).get_id.station_code.toCharArray)' '.' ...
                (jChannels(ii).get_id.site_code.toCharArray)' '.' ...
                (jChannels(ii).get_id.channel_code.toCharArray)'];
        myChannel.name=char(jChannels(ii).name);
        myChannel.startTime=char(jChannels(ii).effective_time.start_time.date_time);
        myChannel.endTime=char(jChannels(ii).effective_time.end_time.date_time);
        myChannel.orientation.azimuth=jChannels(ii).an_orientation.azimuth;
        myChannel.orientation.dip=jChannels(ii).an_orientation.dip;
        myChannel.sampIntv = get_sampIntv(jChannels(ii).sampling_info);
        
        if strcmp(channelOption,'ALL')
            channel_chosen=1;    
        else
            aster_index=find(channelOption=='*');
            code=upper(myChannel.id.channelCode);
            cmp=(channelOption==code);
            cmp(aster_index)=1;
            if all(cmp)
                channel_chosen=1;
            end;
        end;
        if channel_chosen==0 continue;end;
        
        if isempty(availableTime)
            cn=cn+1;
            all_chanNum=all_chanNum+1;
            stationChannel(cn)=myChannel;
            fprintf('  %s',myChannel.id.channelCode);
        else
            if timecmp(availableTime,myChannel.startTime)>=0 & timecmp(availableTime,myChannel.endTime)<0
                cn=cn+1;
                all_chanNum=all_chanNum+1;
                stationChannel(cn)=myChannel;
                fprintf('  %s',myChannel.id.channelCode);
            end;
        end;
    end; %for ii=1:length(jChannels)
    
    if cn > 0
        myNetwork.networkStation(kk).channel=stationChannel;
    end;
    fprintf('   %d channels found\n',cn);
end;

fprintf('%d channels found \n',all_chanNum);
disp('  All done !');
