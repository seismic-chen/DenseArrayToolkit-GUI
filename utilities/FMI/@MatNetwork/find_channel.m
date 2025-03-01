function myNetwork=find_channel(myNetwork,varargin)

% FIND_CHANNEL find channels for network specified by myNetwork.currentNetwork
%
% Usage:
%   myNetwork=find_channel(myNetwork,'criteria',value, .... )
%
% The critieria could be one or combination of fellowing:
%   Criteria Name               Value /[default]
%     Area               [min_lat, max_lat, min_long, max_long] / [-90 90 -180 -180]
%     Sampling:          [min_sampling_rate, max_sampling_rate] / [0.1 1000]
%     Orietation:        [azimuth,dip,angular_distance]
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

if isempty(myNetwork.networkAccess)
    error('Networks not retrieved yet, please run RETRIEVE_NETWORK first');
    return;
end;

if isempty(myNetwork.currentNetwork)
    error('Current_Network is not set. Please run SET_CURRENT_NETWORK first');
    return;
end;

if mod(length(varargin),2)==1
    error('Input arguments not formatted correctly');
end;

%default value
area=[-90 90 -180 180];
sampling=[0.1 1000];
orietation=[0,90,10];

%noCriteria=0;
%if length(varargin)==0
%    noCriteria=1;
%end;

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

netAccess=myNetwork.networkAccess(myNetwork.currentNetwork);
jChannels=netAccess.locate_channels(boxArea,samplingRange,orietationRange);

if length(jChannels)==0
    disp('  No channel found!');
    return;
end;
myNetwork.networkStation=group_by_station(jChannels);

return;

% subroutine group_by_station
function station=group_by_station(jChannels,networkIndex)

station=[];
lastCode=[];
sn=0; % index of station
cn=0; % index of channel
for ii=1:length(jChannels)
    stationCode=char(jChannels(ii).my_site.my_station.get_id.station_code)
    if lastCode~=stationCode
        sn=sn+1;
        lastCode=code;
        station(sn).code=stationCode;
        station(sn).name=(jChannels(ii).my_site.my_station.name.toCharArray);
        station(sn).operator=(jChannels(ii).my_site.my_station.operator.toCharArray);
        station(sn).description=(jChannels(ii).my_site.my_station.description.toCharArray);
        station(sn).comment=(jChannels(ii).my_site.my_station.comment.toCharArray);
        station(sn).startTime=(jChannels(ii).my_site.my_station.effective_time.start_time.time_date.toCharArray)';
        station(sn).endTime=(jChannels(ii).my_site.my_station.effective_time.end_time.time_date.toCharArray)';
        myLocation.latitude=jChannels(ii).my_site.my_station.my_location.latitude;
        myLocation.longitude=jChannels(ii).my_site.my_station.my_location.longitude;
        myLocation.elevation=jChannels(ii).my_site.my_station.my_location.elevation.value;
        station(sn).location=myLocation;
        station(sn).networkIndex=myNetwork.currentNetwork;
        cn=0;
    end;
    myChannel.name=(jChannels(ii).name.toCharArray)';
    myChannel.code=(jChannels(ii).get_id.channel_code.toCharArray);
    myChannel.sampleIntv = get_sampIntv(myNetwork,jChannels(ii).sampling_info);
    myChannel.siteCode=(jChannels(ii).get_id.site_code.toCharArray);
    myChannel.orietation.azimuth=jChannels(ii).an_orietation.azimuth;
    myChannel.orietation.dip=jChannels(ii).an_orietation.dip;
    myChannel.startTime=(jChannels(ii).effective_time.start_time.time_date.toCharArray)';
    myChannel.endTime=(jChannels(ii).effective_time.end_time.time_date.toCharArray)';
    
    channel_duplicate=0;
    for jj=1:cn
        if all(myChannel.code==station(sn).channel(cn).code)
            if timecmp(myChannel.startTime,station(sn).channel(cn).starTime)==1
                station(sn).channel(jj)=myChannel;
            end;
            channel_duplicate=1;
            break;
        end;
    end;
    
    if channel_duplicate==0;
        cn=cn+1;
        station(sn).channel(cn)=myChannel;
    end;
end;
