function [myNetwork, networkNum] = find_AllNetwork(myNetwork)

% GET_ALLNETWORK find all networks and retrieve network infromation
%
% Usage:
%   [myNetwork, networkNum] = find_AllNetwork(myNetwork)
%
% Results:
%   This method will fill three data members, allNetworAccess,
%   allNetworkInfo, and allNetworkNum upon completion. You may use GET
%   function or structure reference to get these data members.
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

disp('Finding networks ... ');
netAccess=myNetwork.networkFnd.retrieve_all;
if length(netAccess)==0
    disp('! Failed finding network.');
    return;
end;
myNetwork.allNetworkAccess=netAccess;
myNetwork.allNetworkNum=length(myNetwork.allNetworkAccess);
disp(sprintf('Totally %d networks found',myNetwork.allNetworkNum));


disp('Retrieving network information for DMC ...');
networkInfo=[];
for ii=1:myNetwork.allNetworkNum
    fprintf('.');
    if mod(ii,40)==0 fprintf('\n');end;
    attrib = netAccess(ii).get_attributes;
    networkInfo(ii).id.networkCode = (attrib.get_id.network_code.toCharArray)';
    networkInfo(ii).id.beginTime = (attrib.get_id.begin_time.date_time.toCharArray)';
    
%     % the following line is just for a temporary fix to a server bug.
%     networkInfo(ii).id.beginTime = [networkInfo(ii).id.beginTime(1:8) 'T' networkInfo(ii).id.beginTime(9:end)];
   
    networkInfo(ii).code = (attrib.get_id.network_code.toCharArray)';
    networkInfo(ii).name = (attrib.name.toCharArray)';
    networkInfo(ii).owner = (attrib.owner.toCharArray)';
    networkInfo(ii).description = (attrib.description.toCharArray)';
    networkInfo(ii).startTime = (attrib.effective_time.start_time.date_time.toCharArray)';
    networkInfo(ii).endTime = (attrib.effective_time.end_time.date_time.toCharArray)';
%     % the following two lines are just for a temporary fix to a server bug.
%     networkInfo(ii).startTime = [networkInfo(ii).startTime(1:8) 'T' networkInfo(ii).startTime(9:end)];
%     networkInfo(ii).endTime = [networkInfo(ii).endTime(1:8) 'T' networkInfo(ii).endTime(9:end)];
end;
myNetwork.allNetworkInfo=networkInfo;
