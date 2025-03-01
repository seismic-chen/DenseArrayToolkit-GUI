function list_station(myNetwork,option)

% LIST_STATION list properties of station
%
% Usage
%   list_stations(myNetwork,'networkStation') list stations in "networkStation"
%   list_stations(myNetwork,'stationPool') list stations in "stationPool".
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

fprintf('      Name    Operator   StartTime    EndTime     Latitude  Longitude   Channels\n')
for ii=1:length(station) 
    channels='';
    for jj=1:length(station(ii).channel)
        channels=sprintf('%s %s',channels,station(ii).channel(jj).code);
    end;
    startTime=string2time(station(ii).startTime);
    st=sprintf('%4d-%02d-%02d',startTime.year,startTime.month,startTime.day);
    endTime=string2time(station(ii).endTime);
    if endTime.year>2100
        et='  N/A     ';
    else
        et=sprintf('%4d-%2d-%2d',endTime.year,endTime.month,endTime.day);
    end;        
    fprintf('%3d  %-7s    %s     %s   %s %8.3f   %8.3f   %s \n',ii,station(ii).name,station(ii).operator, ...
        st,et,station(ii).location.latitude,station(ii).location.longitude,channels);
end;

    