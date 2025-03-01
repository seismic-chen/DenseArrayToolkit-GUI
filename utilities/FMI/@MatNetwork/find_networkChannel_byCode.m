function myNetwork=find_channel_by_code(myNetwork,netId, stationCode, siteCode, channelCode,availableTime)

% FIND_CHANNEL_BY_CODE find channels by station code, site code, and channel code
%
% Usage:
%   myNetwork=find_channel(myNetwork, netId, stationCode, siteCode, channelCode,availableTime)
%
% Input arguments:
%   All codes must be the exact codes. No wildcats will not be recognized.   
%   If availableTime is specified, only the channels that are effective
%   at availableTime will be selected
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

if nargin<4
    error('At least three input arguments needed');
end;

% input arguments checking ...
if exist('availableTime')~=1
    availableTime = [];
end;

if ~isempty(availableTime) & ~isa(availableTime,'char') & ~isa(availableTime,'numeric')
    error('AvailableTime must be Fissures time string or time structure');
end;

if isa(stationCode,'char')
    stationCode={stationCode};
end;
if isa(channelCode,'char')
    channelCode={channelCode};
end;
if isa(siteCode,'char')
    siteCode={siteCode};
end;
stationCode=upper(stationCode);
channelCode=upper(channelCode);
siteCode=upper(siteCode);

jNetId = convert_id(netId,'net');
netAccess = myNetwork.networkFnd.retrieve_by_id(jNetId);

jChannels = netAccess.retrieve_channels_by_code(stationCode, siteCode, channelCode); 

if ~isempty(availableTime)
    timeSelect = [];
    for ii = 1:length(jChannels)
        if timecmp(availableTime,(jChannels(ii).effective_time.start_time.date_time.toCharArray)')>=0 ...
               & timecmp(availableTime,(jChannels(ii).effective_time.end_time.date_time.toCharArray)')<0
            timeSelect = [timeSelect ii];
        end;
    end;
    jChannels = jChannels(timeSelect);
end;

if length(jChannels)==0
    disp('  No channel found!');
    return;
else
    fprintf('  %d channels are found \n', length(jChannels));
end;

myNetwork.networkStation = group_channel_by_station(jChannels);
myNetwork.networkStationNum = length(myNetwork.networkStation);

