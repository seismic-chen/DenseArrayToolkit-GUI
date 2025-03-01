function myNetwork=find_networkChannel(myNetwork,netId,varargin)

% FIND_CHANNEL find channels for network specified by myNetwork.currentNetwork
%
% Usage:
%   myNetwork=find_channel(myNetwork,netId,'criteria',value, .... )
%
% The critieria could be one or combination of fellowing:
%   Criteria Name               Value /[default]
%     Area               [min_lat, max_lat, min_long, max_long] / [-90 90 -180 -180]
%     Sampling:          [min_sampling_rate, max_sampling_rate] / [0.1 1000]
%     Orietation:        [azimuth,dip,angular_distance]
%     AvailableTime:     Only the channels that are effective at
%                        AvailableTime will be selected
%
% If no criteria given, find all stations.
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

if mod(length(varargin),2)==1
    error('Input arguments not formatted correctly');
end;

jNetId = convert_id(netId,'net');
netAccess = myNetwork.networkFnd.retrieve_by_id(jNetId);

%default value
area=[-90 90 -180 180];
sampling=[0.1 1000];
orietation=[0,90,10];
availableTime = [];

% input parameter check
for ii=1:2:length(varargin)
    switch lower(varargin{ii})
    case 'area'
        area=varargin{ii+1};
        if ~isa(area,'numeric') | (isa(area,'numeric') & length(area(:))~=4)
            error('Area must be a 4-element numeric array');
        end;
    case 'sampling'
        sampling=varargin{ii+1};
        if ~isa(sampling,'numeric') | (isa(sampling,'numeric') & length(sampling(:))~=2)
            error('SamplingRange must be a 2-element numeric array');
        end;
    case 'orietation'
        orietation=varargin{ii+1};
        if ~isa(orietation,'numeric') | (isa(orietation,'numeric') & length(orietation(:))~=4)
            error('Area must be a 4-element numeric array');
        end;
    case 'availabletime'
        availableTime = varargin{ii+1};
        if ~(isa(availableTime, 'char') | isa(availableTime, 'numeric'))
            error('AvailableTime must be a string or FMI time stucture');
        end;
    otherwise
        error(['Unrecognized criteria ' varargin{ii}]);
    end;
end;

boxArea=BoxAreaImpl(area(1),area(2),area(3),area(4));

if sampling(1)>=1
    sampMin=SamplingImpl(sampling(1),TimeInterval(1,UnitImpl.SECOND));
else
    sampMin=SamplingImpl(1,TimeInterval(1/sampling(1),UnitImpl.SECOND));
end;
if sampling(2)>=1
    sampMax=SamplingImpl(sampling(2),TimeInterval(1,UnitImpl.SECOND));
else
    sampMax=SamplingImpl(1,TimeInterval(1/sampling(2),UnitImpl.SECOND));
end;
samplingRange=SamplingRangeImpl(sampMin,sampMax);

orCenter=Orientation(orietation(1),orietation(2));
angularDistance=QuantityImpl(orietation(3), UnitImpl.DEGREE);
orietationRange=OrientationRangeImpl(orCenter,angularDistance);

jChannels=netAccess.locate_channels(boxArea,samplingRange,orietationRange);

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

