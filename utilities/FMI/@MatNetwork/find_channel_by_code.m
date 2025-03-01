function myNetwork=find_channel_by_code(myNetwork,stationCode,channelCode,availableTime)

% FIND_CHANNEL_BY_CODE find channels by station code and channel code
%
% Usage:
%   myNetwork=find_channel(stationCode,channelCode,availableTime)
%
% Input arguments:
%   Station code and channel code must be specified. You can search for one code or
%   multiple codes at one time. Specify code as cell array for multiple code search,
%   and string for single code search. '*' can be used for any character match. 
%   For example, '**Z' is for all Z component channels. 'all' can be used for all 
%   stations or all channels
%   
%   "AvailableTime": If availableTime is specified, search channels that are
%       effective at specified time. If it is not specified, search channels that 
%       are effective now
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

if nargin<3
    error('At least three input arguments needed');
end;

% input arguments checking ...
if exist('availableTime')~=1 | isempty(availableTime)
    currentTime=Calendar.getInstance;
    availableTime=sprintf('%4d-%02d-%02dT00:00:00.0Z',currentTime.get(Calendar.YEAR), ...
        currentTime.get(Calendar.MONTH),currentTime.get(Calendar.DATE));
end;
if ~isa(availableTime,'char') & ~isa(availableTime,'struct')
    error('AvailableTime must be Fissures time string or time structure');
end;

if isa(stationCode,'char')
    stationCode={stationCode};
end;
if isa(channelCode,'char')
    channelCode={channelCode};
end;
    
stationCode=upper(stationCode);
channelCode=upper(channelCode);

if isempty(myNetwork.networkAccess)
    error('Networks not retrieved yet, please run RETRIEVE_NETWORK first');
    return;
end;

if isempty(myNetwork.currentNetwork)
    error('Current_Network is not set. Please run SET_CURRENT_NETWORK first');
    return;
end;

networkAccess=myNetwork.networkAccess(myNetwork.currentNetwork);
iterHolder=ChannelIdIterHolder;
disp('Retrieving all channels, please wait ... ');
chanIds=networkAccess.retrieve_all_channels(5000, iterHolder); 
disp('Retrieving channels done!');

if length(chanIds)==0
    disp('  No channel found!');
    return;
end;

% ---------------- Find channels that satisfy search criterias -----------
disp('Filtering channels ... ');
networkStation=[];
sn=0;
myChannel=[];
cn=0;
all_chanNum=0;
for ii=1:length(chanIds)
    
    % -------- first find channels that match input channel id and station id -------
    sta_match=0;
    chan_match=0;
    this_stationCode=(chanIds(ii).station_code.toCharArray)';
    this_channelCode=(chanIds(ii).channel_code.toCharArray)';
    %fprintf('Station: %4s    Channel: %s \n',this_stationCode, this_channelCode);
    if strcmp(stationCode{1},'ALL')
        sta_match=1;
    else
        for jj=1:length(stationCode)
        if length(stationCode{jj})==length(this_stationCode)
            aster_index=find(stationCode{jj}=='*');
            cmp=(stationCode{jj}==this_stationCode);
            cmp(aster_index)=1;
            if all(cmp)
                sta_match=1;
                break;
            end;
        end; % if
        end; % for jj
    end;
    if ~sta_match
        continue;
    end;
    
    if strcmp(channelCode{1},'ALL')
        chan_match=1;
    else
        for jj=1:length(channelCode)
        if length(channelCode{jj})==length(this_channelCode)
            aster_index=find(channelCode{jj}=='*');
            cmp=(channelCode{jj}==this_channelCode);
            cmp(aster_index)=1;
            if all(cmp)
                chan_match=1;
                break;    
            end;
        end; % if
        end; % for jj
    end;
    if ~chan_match
        continue;
    end;
    
    %% ----- match both station code and channel code -------
    jChannel=networkAccess.retrieve_channel(chanIds(ii));
    myChannel.name=char(jChannel.name);
    myChannel.code=char(jChannel.get_id.channel_code);
    myChannel.startTime=char(jChannel.effective_time.start_time.date_time);
    myChannel.endTime=char(jChannel.effective_time.end_time.date_time);
    myChannel.siteCode=char(jChannel.get_id.site_code);
    myChannel.orientation.azimuth=jChannel.an_orientation.azimuth;
    myChannel.orientation.dip=jChannel.an_orientation.dip;

    myChannel.sampIntv = get_sampIntv(myNetwork, jChannel.sampling_info);

    %suppose the output of retrieve_all_channel is sorted by station code
    if sn==0 | (sn>0 & ~strcmp(networkStation(sn).code,this_stationCode))
        % ------------ New Station -------------------
        cn=0;
        jStation=jChannel.my_site.my_station;
        MYC_ST = myChannel.startTime(1:end); ST = MYC_ST;%[yr mth dy] = ymd ( str2num(MYC_ST(1:7)) ); ST = [yr mth dy str2num(MYC_ST(9:10)) str2num(MYC_ST(11:12)) str2num(MYC_ST(13:end))];
        MYC_ET = myChannel.endTime(1:end);   ET = MYC_ET;%[yr mth dy] = ymd ( str2num(MYC_ET(1:7)) ); ET = [yr mth dy str2num(MYC_ET(9:10)) str2num(MYC_ET(11:12)) str2num(MYC_ET(13:end))];
        
        if timecmp(availableTime,time2string(ST))>=0 ...
                & timecmp(availableTime,time2string(ET))<0

            cn=cn+1;
            sn=sn+1;
            all_chanNum=all_chanNum+1;
            %networkStation(sn).code=this_stationCode;
            networkStation(sn).code=char(jStation.get_id.station_code);
            networkStation(sn).name=char(jStation.name);
            networkStation(sn).operator=char(jStation.operator);
            networkStation(sn).description=char(jStation.description);
            networkStation(sn).comment=char(jStation.comment);
            networkStation(sn).startTime=char(jStation.effective_time.start_time.date_time);
            networkStation(sn).endTime=char(jStation.effective_time.end_time.date_time);
            myLocation.latitude=jStation.my_location.latitude;
            myLocation.longitude=jStation.my_location.longitude;
            myLocation.elevation=jStation.my_location.elevation.value;
            networkStation(sn).location=myLocation;
            networkStation(sn).networkCode=myNetwork.networkInfo(myNetwork.currentNetwork).code;
            networkStation(sn).networkIndex=myNetwork.currentNetwork;
            networkStation(sn).channel(cn)=myChannel;
            fprintf('  %-14s    start:%s\n',myChannel.name,myChannel.startTime);
            
        end;
    else
        % -------------- Station already exists -----------------
        MYC_ST = myChannel.startTime(1:end); ST = MYC_ST;
        MYC_ET = myChannel.endTime(1:end);   ET = MYC_ET;
        
        if timecmp(availableTime,time2string(ST))>=0 ...
                & timecmp(availableTime,time2string(ET))<0
        
            cn=cn+1;
            all_chanNum=all_chanNum+1;
            networkStation(sn).channel(cn)=myChannel;
            fprintf('  %-14s    start:%s\n',myChannel.name,myChannel.startTime);
            
        end;
    end;  % if else      
        
end;  %for ii=1:length(chanIds)

myNetwork.networkStation=networkStation;
myNetwork.networkStationNum=length(networkStation);
fprintf('%d stations and %d channels found \n',length(networkStation),all_chanNum);
disp('  All done !');

myNetwork.findStationDone=0;
