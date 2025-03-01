function myNetwork=find_station_channle(myNetwork,channelOption,availableTime)

% FIND_STATION_CHANNEL find channels for stations in field networkStation
%
% Usage: myNetwork=find_station_channel(myNetwork,channel_option,availableTime)
%
% Input arguments:
%   "channel_option" must be string of length 3 to specify components. 
%   Asterisk could be used for substitution of any charactor
%   Example:
%       'BH*':  all components of broad band
%       'BHZ':  z component of broad band
%       '**Z':  z component of all period
%       'ALL':  all channels
%   "AvailableTime": If availableTime is specified, the funciton search channels that are
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

if myNetwork.findStationDone~=1
    disp('Please run FIND_STATION first');
    return;
end;

% input arguments checking ...
if exist('availableTime')~=1 | isempty(availableTime)
    currentTime=Calendar.getInstance;
    availableTime=sprintf('%4d-%02d-%02dT00:00:00.0Z',currentTime.get(Calendar.YEAR), ...
        currentTime.get(Calendar.MONTH),currentTime.get(Calendar.DATE));
end;

if nargin<2
    error('At least 2 input arguments needed');
end;
if ~isa(availableTime,'char') & ~isa(availableTime,'struct')
    error('AvailableTime must be Fissures time string or time structure');
end;
if length(channelOption)~=3 | ~isa(channelOption,'char')
    error('channelOption must be 3 charactors string');
end;
channelOption=upper(channelOption);

% begin query
all_chanNum=0;
currentNetwork=myNetwork.currentNetwork;
netAccess=myNetwork.networkAccess(myNetwork.currentNetwork);
networkId=NetworkId(myNetwork.networkInfo(currentNetwork).code, ...
        edu.iris.Fissures.Time(myNetwork.networkInfo(currentNetwork).startTime,0));
for kk=1:myNetwork.networkStationNum
    fprintf('Retrieving channels for station %4s ... \n',myNetwork.networkStation(kk).code);
    stationId=StationId(networkId,myNetwork.networkStation(kk).code,...
            edu.iris.Fissures.Time(myNetwork.networkStation(kk).startTime,0));
    jChannels=netAccess.retrieve_for_station(stationId);
    cn=0;  % index of channel
    stationChannel=struct([]);
    for ii=1:length(jChannels)
        channel_chosen=0;
        myChannel.name=char(jChannels(ii).name);
        myChannel.code=char(jChannels(ii).get_id.channel_code);
        myChannel.startTime=char(jChannels(ii).effective_time.start_time.date_time);
        myChannel.endTime=char(jChannels(ii).effective_time.end_time.date_time);
        myChannel.siteCode=char(jChannels(ii).get_id.site_code);
        myChannel.orientation.azimuth=jChannels(ii).an_orientation.azimuth;
        myChannel.orientation.dip=jChannels(ii).an_orientation.dip;
%         if jChannels(ii).sampling_info.numPoints>0
%             myChannel.sampIntv=jChannels(ii).sampling_info.interval.value/jChannels(ii).sampling_info.numPoints;
%         else
%             myChannel.sampIntv=nan;
%             disp('Sampling rate is not properly set');
%         end;
        myChannel.sampIntv = get_sampIntv(myNetwork, jChannels(ii).sampling_info);
        
        if strcmp(channelOption,'ALL')
            channel_chosen=1;    
        else
            aster_index=find(channelOption=='*');
            code=upper(myChannel.code);
            cmp=(channelOption==code);
            cmp(aster_index)=1;
            if all(cmp)
                channel_chosen=1;
            end;
        end;
        if channel_chosen==0 continue;end;
        
        % check effective time
        %fprintf('%-4s  %s   Start: %s   end: %s \n', ...
        %        myNetwork.networkStation(kk).code,myChannel.code,myChannel.startTime,myChannel.endTime);
        MYC_ST = myChannel.startTime(1:end-3);ST = [str2num(MYC_ST(1:4)) str2num(MYC_ST(5:6)) str2num(MYC_ST(7:8)) str2num(MYC_ST(9:10)) str2num(MYC_ST(11:12)) str2num(MYC_ST(13:end))];
        MYC_ET = myChannel.endTime(1:end-3);  ET = [str2num(MYC_ET(1:4)) str2num(MYC_ET(5:6)) str2num(MYC_ET(7:8)) str2num(MYC_ET(9:10)) str2num(MYC_ET(11:12)) str2num(MYC_ET(13:end))];
        if timecmp(availableTime,time2string(ST))>=0 ...
                & timecmp(availableTime,time2string(ET))<0
            cn=cn+1;
            all_chanNum=all_chanNum+1;
            if cn==1
                stationChannel=myChannel;
            else
                 stationChannel(cn)=myChannel;
            end;
            fprintf('  %s',myChannel.code);
        end;
        
        % -----------------   Old version   ---------------------
        % remove duplicates
        %for jj=1:cn
        %    if strcmp(myChannel.code,stationChannel(cn).code)
        %        channel_chosen=0;
        %        if timecmp(myChannel.startTime,stationChannel(cn).startTime)>=0
        %        % found duplicate, chosen the newest one
        %            stationChannel(cn)=myChannel;
        %        end;
        %        break;
        %    end;
        %end;
        %
        %if channel_chosen==1;
        %    cn=cn+1;
        %    if cn==1
        %        stationChannel=myChannel;
        %    else
        %         stationChannel(cn)=myChannel;
        %    end;
        %    fprintf('%s  ',myChannel.code);
        %end;
        % -------------------------------------------------------
        
    end; %for ii=1:length(jChannels)
    myNetwork.networkStation(kk).channel=stationChannel;
    fprintf('   %d channels found\n',length(stationChannel));
end;

fprintf('%d channels found \n',all_chanNum);
disp('  All done !');
