function list_network(myNetwork)

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

if isempty(myNetwork.networkInfo)
    disp('No network information found');
    return;
end;

fprintf('     Code  StartTime    EndTime      Name \n')
for ii=1:length(myNetwork.networkInfo) 
    netInfo=myNetwork.networkInfo(ii);
    startTime=string2time(netInfo.startTime);
    st=sprintf('%4d-%02d-%02d',startTime.year,startTime.month,startTime.day);
    endTime=string2time(netInfo.endTime);
    if endTime.year>2100
        et='  N/A     ';
    else
        et=sprintf('%4d-%2d-%2d',endTime.year,endTime.month,endTime.day);
    end;        
    fprintf('%3d  %3s   %s   %s   %s\n',ii,netInfo.code,st,et,netInfo.name);
end;

    