function [mySeis, n]=seis_query(mySeis,timeWindow,network,stations,site,channel)

% SEIS_QUERY check available data for given time window and channels
% 
% mySeis=mySeis=seis_query(mySeis,timeWindow,network,stations,site,channel)
%
% Input arguments:
%  timeWindow:  Time window of request filter. It's a structure containing
%               two fields, startTime and endTime, which could be either
%               staring or Time structur
%     network:  network code. Wildcats is acceptable. But only one string 
%               can be specified.
%    stations:  code of stations. It should be a cell array for
%               multiple stations. It could be '*' for all stations in the
%               network.
%        site:  site code (could be empty)
%     channel:  channel code. Wildcard could be used
%               'BH*':  all components of broad band
%               'BHZ':  z component of broad band
%               '**Z':  z component of all period
%          '*', 'ALL':  all channels
%
% Note that POND server contains two catalogs, FARM and SPYDER. The server
% will return data from FARM if it is available, otherwise from SPYDER 
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   July, 2003
%

import edu.iris.Fissures.IfNetwork.*;
import edu.iris.Fissures.network.*;
import edu.iris.Fissures.IfSeismogramDC.*;
import edu.iris.Fissures.model.*;
import edu.iris.Fissures.utility.*;
import java.io.*;
import java.lang.*;

if nargin<4
    error('timeWindow, network, and stations must be specified');
end;

if exist('site')~=1 || isempty(site)
    site = '';
end;
if exist('channel')~=1 || isempty(channel)
    channel = '*';
end;

if isnumeric(timeWindow.startTime)
    timeWindow.startTime = time2string(timeWindow.startTime);
end;
if isnumeric(timeWindow.endTime)
    timeWindow.endTime = time2string(timeWindow.endTime);
end;

if strcmp(class(stations),'char')
    station_n = 1;
    stations = {stations};
elseif strcmp(class(stations),'cell')
    station_n = length(stations);
else
    error('stations much be string for signal station and cell array for multiple stations');
end;

if ~(strcmp(class(network),'char') && strcmp(class(site),'char') && strcmp(class(channel),'char') )
    error('network, site, and channel must be string');
end;

if strcmp(channel,'ALL') 
    channel_option='*';
end;

seisDc = mySeis.seisDc;

clear requestFilter;
for ii = 1:station_n
    sta = stations{ii};
    chanId = ChannelId(NetworkId(network, edu.iris.Fissures.Time('0000',0)), ...
        sta, site,channel, edu.iris.Fissures.Time('0000',0));
    requestFilter(ii) = RequestFilter(chanId, edu.iris.Fissures.Time(timeWindow.startTime,0), ...
                        edu.iris.Fissures.Time(timeWindow.endTime,0) );
end;

disp('  Searching for available seismograms ... ');
availSeisFilt = seisDc.available_data(requestFilter);

% % detach database tag for station code
% availFrom = {};
% for ii = 1:length(availSeisFilt)
%     sta_temp = char(availSeisFilt(ii).channel_id.station_code);
%     jj=findstr('-',sta_temp);
%     if length(jj)==1;
%         station_code = sta_temp(1:jj-1);
%         availFrom{ii} = sta_temp(jj+1:end);
%     else
%         error('  Assertion error'); % assert that there must be exactly one '-' in station/avail_from name
%     end;
%     availSeisFilt(ii).channel_id.station_code = station_code; 
% end;

% mySeis.availFrom = availFrom;

mySeis.availSeisFilt = availSeisFilt;
mySeis.availChanLoc = [];
mySeis.selectedSeis = [];
n = length(mySeis.availSeisFilt);

% reset retrieve iterator
mySeis = init_iterator(mySeis, 0);

if length(mySeis.availSeisFilt) == 0
    fprintf('  No seismogram found. \n');
else
    fprintf('  %d seismograms found. \n', n);
end;

return;
