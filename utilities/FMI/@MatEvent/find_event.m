function myEvent=find_event(myEvent,varargin)

% FIND_EVENT searches events for given criterias
%
% Usage: myEvent=find_event(myEvent,'PropertyName1',PropertyValue1 ... )
%        myEvent must be an object of class MatEvent
%
% Properties could be any one or more of following
%     Property Name       Property Value / (default value)
%       Area                [min_lat, max_lat, min_long, max_long] / [-90 90 -180 -180]
%       Depth               [min_depth,max_depth] / [0 100]
%       TimeRange                {'start_time','end_time'} / (last month)
%       Magnitude           [min_mag,max_mag] / [5 9]
%       SearchTypes         Type name / {'MB'}
%       Catalogs            Catalog name / {'PREF'}
%       Contributors        Contributor name / {'%'}  - Not used currently
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

import edu.iris.Fissures.IfEvent.*;
import edu.iris.Fissures.*;
import edu.iris.Fissures.model.*;
import edu.iris.Fissures.utility.*;
import java.io.*;
import java.lang.*;
import java.util.*;

msg=nargchk(1,15,nargin);
if ~isempty(msg)
    disp(msg);
    return;
end;

if ~isa(myEvent, 'MatEvent')
    error('myEventDC must be type of MatEvent class');
end;

area=[-90 90 -180 180];
depth=[0 100];
currentTime=Calendar.getInstance;
time{1}=sprintf('%4d-%02d-%02dT00:00:00.0Z',currentTime.get(Calendar.YEAR), ...
    currentTime.get(Calendar.MONTH),currentTime.get(Calendar.DATE));
time{2}=sprintf('%4d-%02d-%02dT00:00:00.0Z',currentTime.get(Calendar.YEAR), ...
    currentTime.get(Calendar.MONTH)+1,currentTime.get(Calendar.DATE));
magnitude=[5.0 8.0];
searchTypes={'MB'};
catalogs={'PREF'};
contributors={'%'};

% input parameter check
for ii=1:2:length(varargin)
    switch lower(varargin{ii})
    case 'area'
        area=varargin{ii+1};
        if ~isa(area,'numeric') | (isa(area,'numeric') & length(area(:))~=4)
            error('Area must be a 4-element numeric array');
        end;
    case 'depth'
        depth=varargin{ii+1};
        if ~isa(depth,'numeric') | (isa(depth,'numeric') & length(depth(:))~=2)
            error('Depth must be a 2-element numeric array');
        end;
    case 'timerange'
        time=varargin{ii+1};
        if ~isa(time,'cell') | (isa(time,'cell') & length(time)~=2)
            error('Time must be a 2-element cell array');
        elseif ~(isa(time{1},'char') & isa(time{2},'char'))
            error('Each element of Time must be a string');
        end;
    case 'magnitude'
        magnitude=varargin{ii+1};
        if ~isa(magnitude,'numeric') | (isa(magnitude,'numeric') & length(magnitude(:))~=2)
            error('Magnitude must be a 2-element numeric array');
        end;
    case 'searchtypes'
        searchTypes=varargin{ii+1};
        if ~isa(searchTypes,'cell') 
            error('SearchTypes must be a cell array');
        end;
    case 'catalogs'
        catalogs=varargin{ii+1};
        if ~isa(catalogs,'cell')
            error('catalogs must be a cell array');
        end;
    otherwise
        error(['Unrecognized parameter ' varargin{ii}]);
    end;
end;
    
boxArea=BoxAreaImpl(area(1),area(2),area(3),area(4));
minDepth=QuantityImpl(depth(1),UnitImpl.KILOMETER);
maxDepth=QuantityImpl(depth(2),UnitImpl.KILOMETER);
timeRange=TimeRange(edu.iris.Fissures.Time(time{1},-1), ...
                    edu.iris.Fissures.Time(time{2},-1));
iterHolder=EventSeqIterHolder;
seq_max=10000;
myEvent.eventAccess=myEvent.eventFnd.query_events(boxArea,minDepth,maxDepth,timeRange, ...
            searchTypes,magnitude(1),magnitude(2),catalogs,contributors,seq_max,iterHolder);
myEvent.nEvent=length(myEvent.eventAccess);
myMatEvent.name = {};
myMatEvent.origin = {};
myMatEvent.preferredOrigin = struct([]);
disp(sprintf('  %d events found\n',myEvent.nEvent));
        
