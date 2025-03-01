function list_channel(myNetwork,sta,option)

% LIST_CHANNEL list channel information for specified station
% 
% Usage:
%   list_channel(myNetwork,station_index,'NetworkStation') 
%   list_channel(myNetwork,station_code,'StationPool')
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

switch lower(option)
case 'networkstation'
    station=myNetwork.networkStation;
case 'stationpool'
    station=myNetwork.stationPool;
otherwise
    error('Unregnized option');
end;
if isempty(station)
    disp('No stations found');
    return;
end;

stationFound=0;
thisStation=[];
if ischar(sta)
    for ii=1:length(station)
        if strcmp(station(ii).code,upper(sta))
            thisStation=station(ii);
            stationFound=1;
            break;
        end;
    end;
elseif isa(sta,'double')
    if sta>0 & sta<=length(station)
        thisStation=station(sta);
        stationFound=1;
    end;
else
    error('Incorrect input argument');
end;

if stationFound==0
    disp('No station found');
    return;
end;

if isempty(thisStation)
    disp('No channel found in this station');
end;

channel=thisStation.channel;

fprintf('        Name         StartTime    EndTime     Azimuth     Dip   SampleIntv \n')
for ii=1:length(channel)
    startTime=string2time(channel(ii).startTime);
    st=sprintf('%4d-%02d-%02d',startTime.year,startTime.month,startTime.day);
    endTime=string2time(channel(ii).endTime);
    if endTime.year>2100
        et='  N/A     ';
    else
        et=sprintf('%4d-%2d-%2d',endTime.year,endTime.month,endTime.day);
    end;        
    fprintf('%2d  %-7s    %s   %s   %5.1f    %5.1f    %f \n',ii,channel(ii).name, ...
        st,et,channel(ii).orientation.azimuth,channel(ii).orientation.dip,channel(ii).sampIntv);
end;