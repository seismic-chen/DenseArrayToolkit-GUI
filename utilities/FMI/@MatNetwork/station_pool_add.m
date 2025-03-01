function myNetwork=station_pool_add(myNetwork);

% STATION_POOL_ADD add stations to station pool
%
% Usage:
%   myNetwork=station_pool_add(myNetwork) add stations in field "networkStation" 
%           station pool.
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

if myNetwork.networkStationNum==0
    disp('No station found in station pool');
    return;
end;

n=myNetwork.stationPoolNum;
for ii=1:myNetwork.networkStationNum
    if ~isempty(myNetwork.networkStation(ii).channel)
        n=n+1;
        if n==1
            myNetwork.stationPool=myNetwork.networkStation(ii);
        else
            myNetwork.stationPool(n)=myNetwork.networkStation(ii);
        end;
    end;
end;

myNetwork.stationPoolNum=length(myNetwork.stationPool);
myNetwork.stationAddedToPool=1;
