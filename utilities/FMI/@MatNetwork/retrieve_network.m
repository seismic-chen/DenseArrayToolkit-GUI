function myNetwork=retrieve_network(myNetwork)

% RETRIEVE_NETWORK retrieves all network from DMC
%
% Usage:
%   myNetwork=retrieve_network(myNetwork)
%
% Note:
%   There are two steps. Firstly, only obtain accesses to networks;
%   Secondly, retrieve informations for all networks. Since the second step 
%   is slow, so an information data file "networkInfo.mat" is saved. Next time
%   the function will locate this file before retrieving information from DMC. 
%   If the informaion file exists and there is no new network added since last 
%   time, the funciton will load information into the object. Otherwise the 
%   function still retrieves informations from DMC.
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

if ~isempty(myNetwork.networkAccess)
    disp('Networks are already retrieved');
    return;
end;

disp('Finding networks ... ');
netAccess=myNetwork.networkFnd.retrieve_all;
if length(netAccess)==0
    disp('! Failed finding network.');
    return;
end;
myNetwork.networkAccess=netAccess;
myNetwork.networkNum=length(myNetwork.networkAccess);
disp(sprintf('Totally %d networks found',myNetwork.networkNum));

if exist('networkInfo.mat')==2
    fname=which('networkInfo.mat');
    eval(['load ' fname]);
    if length(networkInfo)==myNetwork.networkNum
        disp(['Network information file loaded: ' fname])
        myNetwork.networkInfo=networkInfo;
    else
        disp(['Network information file seems to be out of date: ' fname]);
        disp('Retrieving network information for DMC');
        disp('Please wait ... ');
        
        networkInfo=[];
        for ii=1:myNetwork.networkNum
            fprintf('.');
            if mod(ii,40)==0 fprintf('\n');end;
            attrib = netAccess(ii).get_attributes;
            networkInfo(ii).code=(attrib.get_id.network_code.toCharArray)';
            networkInfo(ii).name=(attrib.name.toCharArray)';
            networkInfo(ii).owner=(attrib.owner.toCharArray)';
            networkInfo(ii).description=(attrib.description.toCharArray)';
            networkInfo(ii).startTime=(attrib.effective_time.start_time.date_time.toCharArray)';
            networkInfo(ii).endTime=(attrib.effective_time.end_time.date_time.toCharArray)';
        end;
        myNetwork.networkInfo=networkInfo;
        %eval(['save ' fname ' networkInfo']);
    end;
else
        disp('Retrieving network information for DMC');
        disp('Please wait ... ');
        fname=which('open_orb.m');
        fname=strrep(fname,'open_orb.m','networkInfo.mat');
        
        networkInfo=[];
        for ii=1:myNetwork.networkNum
            fprintf('.');
            if mod(ii,40)==0 fprintf('\n');end;
            attrib = netAccess(ii).get_attributes;
            networkInfo(ii).code=(attrib.get_id.network_code.toCharArray)';
            networkInfo(ii).name=(attrib.name.toCharArray)';
            networkInfo(ii).owner=(attrib.owner.toCharArray)';
            networkInfo(ii).description=(attrib.description.toCharArray)';
            networkInfo(ii).startTime=(attrib.effective_time.start_time.date_time.toCharArray)';
            networkInfo(ii).endTime=(attrib.effective_time.end_time.date_time.toCharArray)';
        end;
        myNetwork.networkInfo=networkInfo;
        %eval(['save ' fname ' networkInfo']);
end;
    
