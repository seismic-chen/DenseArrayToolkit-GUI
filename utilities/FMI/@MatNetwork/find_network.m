function myNetwork = find_network(myNetwork, in_str, option );

% FIND_NETWORK find networks by name or code. If all network information
% has been retrieved previously, the function will search allNetworkInfo.
% Other wise the function will search from remote data center.
%
% myNetwork = find_network(myNetwork, in_str, option)
%   Input Arguments:
%      in_str:   input string that may be name or code of the network.
%      option:   'name' or 'code'.
%
% Results:
%   Three data members, networAccess, networkInfo, and networkNum will be
%   set upon completion. You may use GET function or structure reference 
%   to get these data members.
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

if exist('option')~=1 || isempty(option)
    option = 'code';
end;

myNetworkDC.networkInfo=struct([]);
myNetworkDC.networkAccess=[];
myNetworkDC.networkNum=0;
disp('Finding networks ... ');
if myNetwork.allNetworkNum == 0   % all network inform not available
    switch lower(option)
        case 'code'
            netAccess = myNetwork.networkFnd.retrieve_by_code(in_str);
        case 'name'
            netAccess = myNetwork.networkFnd.retrieve_by_name(in_str);
        otherwise
            error('Unrecognized option');
    end;
    if length(netAccess)==0
        disp('No network found !');
        return;
    end;
    myNetwork.networkAccess=netAccess;
    myNetwork.networkNum=length(myNetwork.networkAccess);
    fprintf('Totally %d networks found \n',myNetwork.networkNum);
    % retrieving network infor
    networkInfo=[];
    for ii=1:myNetwork.networkNum
        attrib = netAccess(ii).get_attributes;
        networkInfo(ii).id.networkCode = (attrib.get_id.network_code.toCharArray)';
        networkInfo(ii).id.beginTime = (attrib.get_id.begin_time.date_time.toCharArray)';
    
%         % the following line is just for a temporary fix to a server bug.
%         networkInfo(ii).id.beginTime = [networkInfo(ii).id.beginTime(1:8) 'T' networkInfo(ii).id.beginTime(9:end)];
   
        networkInfo(ii).code = (attrib.get_id.network_code.toCharArray)';
        networkInfo(ii).name = (attrib.name.toCharArray)';
        networkInfo(ii).owner = (attrib.owner.toCharArray)';
        networkInfo(ii).description = (attrib.description.toCharArray)';
        networkInfo(ii).startTime = (attrib.effective_time.start_time.date_time.toCharArray)';
        networkInfo(ii).endTime = (attrib.effective_time.end_time.date_time.toCharArray)';
%         % following two lines are just for a temporary fix to a server bug.
%         networkInfo(ii).startTime = [networkInfo(ii).startTime(1:8) 'T' networkInfo(ii).startTime(9:end)];
%         networkInfo(ii).endTime = [networkInfo(ii).endTime(1:8) 'T' networkInfo(ii).endTime(9:end)];
    end;
    myNetwork.networkInfo=networkInfo;
else
    found_index = [];
    for ii = 1:myNetwork.allNetworkNum
        switch lower(option)
            case 'code'
                if strcmp(in_str, myNetwork.allNetworkInfo(ii).code)
                    found_index = [found_index ii];
                end;
            case 'name'
                if strcmp(in_str, myNetwork.allNetworkInfo(ii).name)
                    found_index = [found_index ii];
                end;
            otherwise
                error('Unrecognized option');
        end;
    end;        
    if length(found_index)==0
        disp('No network found !');
        return;
    end;
    myNetwork.networkNum = length(found_index);
    myNetwork.networkAccess = myNetwork.allNetworkAccess(found_index);
    myNetwork.networkInfo = myNetwork.allNetworkInfo(found_index);
    fprintf('Totally %d networks found \n',myNetwork.networkNum);
end;
