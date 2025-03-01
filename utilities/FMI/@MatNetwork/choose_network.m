function myNetwork=choose_network(myNetwork,varargin);

% CHOOSE_NETWORK choose networks by name, code, start time, or time range
%
% myNetwork=choose_network(myNetwork,'critieria',value, ...)
%
% The criteria(s) could be:
% Criteria         Description
%  Name:        the name of the network, string for single network and cell
%               array for multiple networks
%  Code:        the code of the network, string for single network and cell
%               array for multiple networks
%  StartTime:   the start time of the network. The code is not unique, but
%               the combination of code nad start time is unique.
%  AvailableTime:   to find the networks that are available at the specified time.
%               AvailableTime could be a string or time structure
%
% If no criteria given, choose all networks.
%
% The indices of chosen networks saved in field chosenNetwork
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

if isempty(myNetwork.networkAccess)
    error('Networks not retrieved yet, please run RETRIEVE_NETWORK first');
end;

if length(varargin)==0
    myNetwork.chosenNetwork=1:myNetwork.networkNum;
end;

if mod(length(varargin),2)==1
    error('Input arguments not formatted correctly');
end;

n=0;
myNetwork.chosenNetwork=[];
for kk=1:length(myNetwork.networkAccess)
    misMatch=0;
    for ii=1:2:length(varargin)
        switch lower(varargin{ii})
        case 'name'
            name=varargin{ii+1};
            if ischar(name)
                name={name};
            end;
            nameMatch=0;
            for jj=1:length(name)
                if strcmp(myNetwork.networkInfo(kk).name,name{jj})
                    nameMatch=1;
                    break;
                end;
            end;
            if nameMatch==0
                misMatch=1;
                continue
            end;
        case 'code'
            code=varargin{ii+1};
            if ischar(code)
                code={code};
            end;
            codeMatch=0;
            for jj=1:length(code)
                if strcmp(myNetwork.networkInfo(kk).code,code{jj})
                    codeMatch=1;
                    break;
                end;
            end;
            if codeMatch==0
                misMatch=1;
                continue
            end;
        case 'starttime'
            if cmp(myNetwork.networkInfo(kk).startTime,varargin{ii+1})~=0
                misMatch=1;
                continue
            end;
        case 'availabletime'
            s=strrep(varargin{ii+1},'-','');
            availableTime=str2num(s(1:8));
            s=strrep(myNetwork.networkInfo(kk).startTime,'-','');
            startTime=str2num(s(1:8));
            s=strrep(myNetwork.networkInfo(kk).endTime,'-','');
            endTime=str2num(s(1:8));
            
            if availableTime<startTime | availableTime>endTime
                misMatch=1;
                continue
            end;
        end;
    end;
    
    if misMatch==0  % match
        n=n+1;
        myNetwork.chosenNetwork(n)=kk;
    end;
end;
            
if n==0
    disp('No network found');
else
    disp(sprintf('  %d networks found',n));
end;

