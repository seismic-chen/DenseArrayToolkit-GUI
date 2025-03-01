function out=subsref(a,s)
% SUBSREF class MatNetwork subscripting expression interpreter
%
% Usage:
%   myNetworkDC.networkNum
%              .networkInfo
%              .networkStationNum
%              .networkStation
%              .stationPoolNum
%              .stationPool
%              .chosennetwork
%              .currentNetwork
%  where fields 'network', 'networkStation', and 'stationPool' could be followed
%     by index, e.g. myNetworkDC.network(2)
%
%  myNetworkDC(i) is allowed and equivalent to myNetworkDC.networkInfo(i)
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%
    
out=[];

switch s(1).type
case '.'
    switch lower(s(1).subs)
    case 'networknum'
        out=a.networkNum;
    case 'currentnetwork'
        out=a.currentNetwork;
        if isempty(a.currentNetwork)
            disp('warning: currentNetwork not set yet');
        end
    case 'chosennetwork'
        out=a.chosenNetwork;
        if isempty(a.chosenNetwork)
            warning('No chosen network found');
        end
    case 'networkstationnum'
        out=a.networkStationNum;
    case 'stationpoolnum'
        out=a.stationPoolNum;
    case 'networkinfo'
        if length(s)==1
            out=a.networkInfo;
        else
            if s(2).type~='()'
                disp('Unsupported subscripting expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension subscripting not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            if a.networkNum==0
                disp('No network found');
                return;
            end;
            if ind<1 | ind>a.networkNum
                disp('Index out of limit');
                disp_usage;
                return;
            end;
            out=a.networkInfo(ind);
        end;
    case    'networkstation'
        if length(s)==1
            out=a.networkStation;
        else
            if s(2).type~='()'
                disp('Unsupported subscripting expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension subscripting not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            if a.networkStationNum==0
                disp('No network station found');
                return;
            end;
            if ind<1 | ind>a.networkStationNum
                disp('Index out of limit');
                disp_usage;
                return;
            end;
            out=a.networkStation(ind);
        end;
    case  'stationpool'
        if length(s)==1
            out=a.stationPool;
        else
            if s(2).type~='()'
                disp('Unsupported subscripting expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension subscripting not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            if a.stationPoolNum==0
                disp('No network station found');
                return;
            end;
            if ind<1 | ind>a.stationPoolNum
                disp('Index out of limit');
                disp_usage;
                return;
            end;
            out=a.stationPool(ind);
        end;
    otherwise
        disp('Unsupported subscripting expression');
        disp_usage;
        return;
    end;
case '()'
    if length(s(1).subs)~=1
        disp('Multiple dimension subscripting not supported');
        disp_usage;
        return;
    end;
    ind=s(1).subs{:};
    if a.networkNum==0
        disp('No network found');
        return;
    end;
    if ind<1 | ind>a.networkNum
        disp('Index out of limit');
        disp_usage;
        return;
    end;
    out=a.networkInfo(ind);
otherwise
    disp('Unsupported subscripting expression');
    disp_usage;
end;

function disp_usage
    disp(' ');
    disp('Class MatNetwork subscripting usage:');
    disp('   myNetworkDC.networkNum');
    disp('              .networkInfo');
    disp('              .networkStationNum');
    disp('              .networkStation')
    disp('              .stationPoolNum')
    disp('              .stationPool')
    disp('              .currentNetwork')
    disp('  where fields "network", "networkStation", and "stationPool" could be followed')
    disp('  by index, e.g. myNetworkDC.network(2)')
    disp(' ');
    disp('  myNetworkDC(i) is allowed and equivalent to myNetworkDC.network(i)')
    disp(' ');

