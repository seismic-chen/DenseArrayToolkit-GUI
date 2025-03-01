function [mySeis,nSeisRetrieved] = seis_retrieve(mySeis, option)

% SEIS_RETRIEVE retreive seismograms from IRIS data center
%
% [mySeis, nSeisRetrieved] = seis_retrieve(mySeis,option) 
% Input arguments:
%    option: 'available' for retrieving all available seismograms, or
%            'selected' for retrieving selected seismograms only. Default
%            is 'available'
% Ouput arguments:
%    nSeisRetrieved:  number of retrieved seismograms. Note that might be
%    different from the number of mySeis.iteraotr.nPerRetrieve, because one
%    request filter may return multiple seimograms if the request windows
%    contains multiple events.
%
% Note: An iterator is implemented to prevent memory overflow when you 
% retrieve a big number of seismograms. When iterator is enabled, it will
% allows users to break the whole retrieving job into small pieces. Users 
% can retrieve a small number of seimograms at one time, which can be 
% specified by method init_iterator. You must call init_iterator first
% before start a new multiple-retrieving job.
%
% However, be aware of that when you do next retrieving, the last seismograms
% saved in the current MatSeismogram object will be cleared. So make sure to 
% save them onto disks before the next retrieving. Also, if you need to 
% retrieve channel information, you have to retrieve it everytime after you
% the seismogram retrieving.
%
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   July, 2003
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also seis_query, seis_select, enable_iterator, init_iterator, disable_iterator
%

import edu.iris.Fissures.IfSeismogramDC.*;
import edu.iris.Fissures.IfNetwork.*;


if exist('option')~=1 || isempty(option)
    option = 'available';
end;

option = lower(option);
if strcmp(option,'available')==0 && strcmp(option,'selected')==0
    error('  Unrecognized option');
end;
if isempty(mySeis.availSeisFilt)
    error('  Please run seis_query first');
end;
if strcmp(option,'selected') & isempty(mySeis.selectedSeis)
    error('  Please run seis_select first');
end;

if strcmp(option,'selected');
    allSeis = mySeis.selectedSeis;
else
    allSeis = 1:length(mySeis.availSeisFilt);
end;

nSeisRetrieved = 0;

if mySeis.iterator.enabled == 0
    currentSelected = allSeis;
else  % using iterator
    if mySeis.iterator.nPerRetrieve == 0
        error('  Please initialize the iterator first');
    end;
    % first time retrieving
    if mySeis.iterator.nIter==0 & mySeis.iterator.nTotal == 0
        mySeis.iterator.nTotal = length(allSeis);
        mySeis.iterator.nRemained = length(allSeis);
        mySeis.iterator.nRetrieved = 0;
    end;
    % check if retrieving has finished
    if mySeis.iterator.nRemained ==0
        error('  No more seismograms to be retrieved');
        return;
    end
    if mySeis.iterator.nPerRetrieve >= mySeis.iterator.nRemained
        currentSelected = allSeis(end-mySeis.iterator.nRemained+1:end);
        nIterRetrieved = mySeis.iterator.nRetrieved + mySeis.iterator.nRemained;
        nIterRemained = 0;
    else
        currentSelected = allSeis(end-mySeis.iterator.nRemained + (1:mySeis.iterator.nPerRetrieve));
        nIterRetrieved = mySeis.iterator.nRetrieved + mySeis.iterator.nPerRetrieve;
        nIterRemained = mySeis.iterator.nRemained - mySeis.iterator.nPerRetrieve;
    end;
end
% create java array for RequestFilter
rf = javaArray('edu.iris.Fissures.IfSeismogramDC.RequestFilter',length(currentSelected));
rf(1:length(currentSelected)) = mySeis.availSeisFilt(currentSelected);

% retrieve seismograms
if length(mySeis.availSeisFilt) > 0
    disp('Retrieving seismograms ...  ');
    seismograms = mySeis.seisDc.retrieve_seismograms(rf);
else
    error('  AvailSeisFilt is empty. Run seis_query first');
end;
% get number of seismograms
mySeis.nSeismogram = length(seismograms);

% initializing structure for seismogram information
mySeis.seismogram = [];
mySeis.seismogramAttr = struct([]);
mySeis.seisChanId = [];
mySeis.seisChanInfo = [];

% extracting seismogram information from java object
seisAttr = [];
seisChanId = javaArray('edu.iris.Fissures.IfNetwork.ChannelId',1);
for ii=1:mySeis.nSeismogram
    if seismograms(ii).can_convert_to_double
        s = double(seismograms(ii).get_as_doubles);
    elseif seismograms(ii).can_convert_to_long
        s = double(seismograms(ii).get_as_longs);
    elseif seismograms(ii).can_convert_to_float
        s = double(seismograms(ii).get_as_floats);
    elseif seismograms(ii).can_convert_to_short
        s = double(seismograms(ii).get_as_shorts);
    else
        error('  Cannot convert ');
    end;
    mySeis.seismogram{ii} = s;
    seisChanId(ii) = seismograms(ii).channel_id;
    
    seisAttr(ii).beginTime = char(seismograms(ii).begin_time.date_time);
    seisAttr(ii).nPoint = seismograms(ii).num_points;
    seisAttr(ii).network = char(seismograms(ii).channel_id.network_id.network_code);
    seisAttr(ii).station = char(seismograms(ii).channel_id.station_code);
    seisAttr(ii).site = char(seismograms(ii).channel_id.site_code);
    seisAttr(ii).channel = char(seismograms(ii).channel_id.channel_code);
    seisAttr(ii).sampleIntv = seismograms(ii).sampling_info.interval.value / ...
                          seismograms(ii).sampling_info.numPoints;
    seisAttr(ii).sampleIntvUnit = char(seismograms(ii).sampling_info.interval.the_units.toString);
    seisAttr(ii).unit = char(seismograms(ii).y_unit.toString);
    if isfield(seismograms(ii),'properties')
        seisAttr(ii).from = char(seismograms(ii).properties(1).value);
        seisAttr(ii).qualityFlag = char(seismograms(ii).properties(2).value);
        seisAttr(ii).event = char(seismograms(ii).properties(3).value);
        seisAttr(ii).seedFrac = char(seismograms(ii).properties(4).value);
    end
end;

if mySeis.nSeismogram>0
    mySeis.seismogramAttr = seisAttr;
    mySeis.seisChanId = seisChanId;
else
    mySeis.seismogramAttr = struct([]);
    mySeis.seisChanId = [];
end;

nSeisRetrieved = mySeis.nSeismogram;

% update retrieve iterator
if mySeis.iterator.enabled == 1
    mySeis.iterator.nRemained = nIterRemained;
    mySeis.iterator.nRetrieved = nIterRetrieved;
    mySeis.iterator.nIter = mySeis.iterator.nIter + 1;
end;

fprintf('  %d seismograms retrieved. \n',mySeis.nSeismogram);

if mySeis.iterator.enabled == 0
    fprintf(' Done! \n');
end
if mySeis.iterator.enabled == 1 & mySeis.iterator.nRemained == 0
    fprintf(' Done! \n');
end
    

