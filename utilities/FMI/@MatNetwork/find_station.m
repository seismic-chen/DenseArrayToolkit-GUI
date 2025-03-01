function myNetwork=find_station(myNetwork,varargin)

% FIND_STATION find stations for network specified by MatNetwork.currentNetwork
%
% The criteria could be:
% Criteria:         Description
%    Code:          String of code for one stations, or cell array of codes of multiple stations.  
%    Area:          [min_lat,max_lat,min_long,max_long] / [all area]
%    AvailableTime: String of FISSURES time or time structure, used to find 
%                       if the networks is effective for specified time
% 
% Found stations are stored in field "networkStation". Field "findStationDone"
% is set to be 1. 
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

if isempty(myNetwork.currentNetwork)
    error('Current_Network is not set. Please run SET_CURRENT_NETWORK first');
    return;
end;

if mod(length(varargin),2)==1
    error('Input arguments not formatted correctly');
end;

currentNetwork=myNetwork.currentNetwork;
fprintf('Retrieving all stations for network %s ...', ...
     myNetwork.networkInfo(currentNetwork).code); 
allStation=myNetwork.networkAccess(currentNetwork).retrieve_stations;
fprintf('  Done \n');

n=0;
station=[];
for ii=1:length(allStation)
    mismatch=0;
    for jj=1:2:length(varargin)
        switch lower(varargin{jj})
        case 'code'
            code=varargin{jj+1};
            if ~isa(code,'cell') & ~isa(code,'char')
                error('"Code" must be a string or cell array');
            end;
            if isa(code,'char')
                code={code};
            end;

            matchFlag=0;
            for kk=1:length(code)
                if strcmp(allStation(ii).get_id.station_code,code{kk})
                    matchFlag=1;
                    break
                end
            end;
            if matchFlag
                mismatch=0;
            else
                mismatch=1;
                break
            end;
        case 'area'
            staLat=allStation(ii).my_location.latitude;
            staLong=allStation(ii).my_location.longitude;
            area=varargin{jj+1};
            if ~isa(area,'numeric') | (isa(area,'numeric') & length(area)~=4)
                error('Value of AREA is not properly fomatted');
            end;
            if staLat<area(1) | staLat>area(2) | staLong<area(3) | staLong>area(4)
                mismatch=1;
                break
            end;
        case 'availabletime'
            availableTime=varargin{jj+1};
            startTime=(allStation(ii).effective_time.start_time.date_time.toCharArray)';
            endTime=(allStation(ii).effective_time.end_time.date_time.toCharArray)';
            
            if ~(timecmp(availableTime,startTime)>=0 & timecmp(availableTime,endTime)<0)
                mismatch=1;
                break;
            end;
        otherwise
            error(sprintf('Unrecognized option %s',varargin{ii}));
        end;
    end; % for jj
        
    if mismatch==0
        n=n+1;
        station(n).code=(allStation(ii).get_id.station_code.toCharArray)';
        station(n).name=(allStation(ii).name.toCharArray)';
        station(n).networkCode=myNetwork.networkInfo(currentNetwork).code;
        station(n).operator=(allStation(ii).operator.toCharArray)';
        station(n).description=(allStation(ii).description.toCharArray)';
        station(n).comment=(allStation(ii).comment.toCharArray)';
        station(n).startTime=(allStation(ii).effective_time.start_time.date_time.toCharArray)';
        station(n).endTime=(allStation(ii).effective_time.end_time.date_time.toCharArray)';
        myLocation.latitude=allStation(ii).my_location.latitude;
        myLocation.longitude=allStation(ii).my_location.longitude;
        myLocation.elevation=allStation(ii).my_location.elevation.value;
        station(n).location=myLocation;
        station(n).networkIndex=myNetwork.currentNetwork;
        station(n).channel=struct([]);
    end;
end;  % for ii

myNetwork.networkStation=station;
myNetwork.networkStationNum=length(station);
myNetwork.findStationDone=1;

fprintf('%d stations found for %s network\n',myNetwork.networkStationNum, ...
        myNetwork.networkInfo(myNetwork.currentNetwork).code);

