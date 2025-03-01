function [myNetwork nChannel] = find_channels(myNetwork,searchMode,varargin)

% FIND_CHANNELS find channels regardless networks and stations
%
% Usage:
%   myNetwork = find_channels(myNetwork,maxChannel,'criteria',value, .... )
%
% The critieria could be one or combination of fellowing:
%   Criteria Name               Value /[default]
%     Area               [min_lat, max_lat, min_long, max_long] / [-90 90 -180 -180]
%     Sampling:          [min_sampling_rate, max_sampling_rate] / [0.1 1000]
%     Orietation:        [azimuth,dip,angular_distance]
%     AvailableTime:     Only the channels that are effective at
%                        AvailableTime will be selected
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

if exist('searchMode')~=1 | isempty(searchMode)
    searchMode = 'query';
end;

if mod(length(varargin),2)==1
    error('Input arguments not formatted correctly');
end;

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

maxChannels = 10;
iterHolder = ChannelIdIterHolder;
howManyRemain = 0;
channels = struct([]);
nChannel = 0;
if strcmp(lower(searchMode), 'query')
    jChanId = myNetwork.networkExp.locate_channels(boxArea,samplingRange,orietationRange,1,iterHolder);
    nChannel = iterHolder.value.how_many_remain + 1;
    fprintf('%d channels found. Channels are not retrieved.\n', nChannel);
else
    jChanId = myNetwork.networkExp.locate_channels(boxArea,samplingRange,orietationRange,1,iterHolder);
    if isempty(jChanId)
        return;
    end;
    allChannel = iterHolder.value.how_many_remain + 1;
    if allChannel > 180
        fprintf('%d channels found. It may takes more than %d minutes to retrieve them all.\n', allChannel, ceil(allChannel/2/60));
        answer = input('Do you want to continue (y/n)? ', 's');
        if strcmp(lower(answer), 'n') | strcmp(lower(answer), 'no')
            return;
        end;
    else
        fprintf('%d channels found. It may take more than %d seconds to retrieve them all.\n', allChannel, ceil(allChannel/2));
    end;
    disp(['Retrieving (1 dot denotes ' num2str(maxChannels) ' channels) ...']);
    jChanId = myNetwork.networkExp.locate_channels(boxArea,samplingRange,orietationRange,maxChannels,iterHolder);
    nChannel = length(jChanId);
    iter = ChannelIdIterHelper.narrow(iterHolder.value);
    seqHolder = ChannelIdSeqHolder(jChanId);
    howManyRemain = iterHolder.value.how_many_remain;
    fprintf('.');
    moreToRetrieve = 1;
 
    fullCodeAll = {};
    while (moreToRetrieve)
        moreToRetrieve = iter.next_n(maxChannels, seqHolder);
        howManyRemain = iter.how_many_remain;
        this_jChanId = seqHolder.value;
        jChanId = [jChanId this_jChanId];
        nChannel = nChannel + length(this_jChanId);
        fprintf('.');
        if mod(nChannel/maxChannels, 40) == 0
            fprintf(' %d%%  %d remain \n', round(nChannel/allChannel*100), howManyRemain);
        end;
        
        % for debugging
        if howManyRemain < maxChannels*2 & moreToRetrieve
        %if nChannel>1000
            for kk = 1:length(this_jChanId)
                fullCode = [(this_jChanId(kk).network_id.network_code.toCharArray)' '.' ...
                        (this_jChanId(kk).station_code.toCharArray)' '.' ...
                        (this_jChanId(kk).site_code.toCharArray)' '.' ...
                        (this_jChanId(kk).channel_code.toCharArray)' '    ' ...
                        (this_jChanId(kk).begin_time.date_time.toCharArray)'];
                %fullCodeAll{nChannel-maxChannels + kk} = fullCode; 
                fprintf('%s \n',fullCode);
            end;
        end;
    end;
    fprintf('  Done!\n');
end;

