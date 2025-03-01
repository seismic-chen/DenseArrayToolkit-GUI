function out=subsref(myEvent,s)
%
% SUBSREF Method for MatEvent object
%
% Usage:
%   myEvent.nEvent                - total number of events
%          .name(i)               - name of (ith) event 
%          .origin(i)             - origins of (ith) event
%          .preferredOrigin{i}    - preferred origin of (ith) event
%          .preferredOrigin{i}.magnitude
%          .preferredOrigin{i}.location
%          .preferredOrigin{i}.originTime
%          .preferredOrigin{i}.catalog
%          .preferredOrigin{i}.contributor
%  
% The index is optional. It may be a number or a vector.
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

out=[];

if strcmp(s(1).type, '.')~=1
    disp('Unsupported subscripted reference');
    disp_usage;
    return;
end;

switch lower(s(1).subs)
    case 'nevent'
        out=myEvent.nEvent;
    case 'name'
        if isempty(myEvent.name)
            beep;
            disp('Warning: Event informations are not retrieved yet!')
        end;
        if length(s)==1
            out=myEvent.name;
        else
            if s(2).type~='.'
                disp('Unsupported subscripted reference');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension subscripting not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            for kk = 1:length(ind)
                out{kk} = myEvent.name{ind(kk)};
            end;
        end;
    case 'preferredorigin'
        if isempty(myEvent.preferredOrigin)
            beep;
            disp('Warning: Event informations are not retrieved yet!')
        end;
        if length(s)==1
            out=myEvent.preferredOrigin;
        elseif length(s) == 2
            if strcmp(s(2).type,'()') | strcmp(s(2).type,'{}')
                if length(s(2).subs)~=1
                    disp('Multiple dimension subscripting not supported');
                    disp_usage;
                    return;
                end;
                ind=s(2).subs{:};
                for kk = 1:length(ind)
                    out{kk} = myEvent.preferredOrigin(ind(kk));
                end;
                if length(ind)==1
                    out = out{1};
                end;
            else
                switch lower(s(2).subs)
                    case 'magnitude'
                        for kk = 1:myEvent.nEvent
                            out{kk} = myEvent.preferredOrigin(kk).magnitude;
                        end;
                    case 'location'
                        clear out;
                        for kk = 1:myEvent.nEvent
                            out(kk) = myEvent.preferredOrigin(kk).location;
                        end;
                    case 'origintime'
                        for kk = 1:myEvent.nEvent
                            out{kk} = myEvent.preferredOrigin(kk).originTime;
                        end;
                   case 'catalog'
                        for kk = 1:myEvent.nEvent
                            out{kk} = myEvent.preferredOrigin(kk).catalog;
                        end;
                    case 'contributor'
                        for kk = 1:myEvent.nEvent
                            out{kk} = myEvent.preferredOrigin(kk).contributor;
                        end;
                    otherwise
                        disp('Unsupported subscripted reference');
                        disp_usage;
                        return;
                end;
            end;
        elseif length(s) == 3
            if ~(strcmp(s(2).type,'()') | strcmp(s(2).type,'{}')) | s(3).type~='.' 
                disp('Unsupported subscripted reference');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension subscripting not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            switch lower(s(3).subs)
                case 'magnitude'
                    for kk = 1:length(ind)
                        out{kk} = myEvent.preferredOrigin(ind(kk)).magnitude;
                    end;
                case 'location'
                    clear out
                    for kk = 1:length(ind)
                        out(kk) = myEvent.preferredOrigin(ind(kk)).location;
                    end;
                case 'origintime'
                    for kk = 1:length(ind)
                        out{kk} = myEvent.preferredOrigin(ind(kk)).originTime;
                    end;
                case 'catalog'
                    for kk = 1:length(ind)
                        out{kk} = myEvent.preferredOrigin(ind(kk)).catalog;
                    end;
                case 'contributor'
                    for kk = 1:length(ind)
                        out{kk} = myEvent.preferredOrigin(ind(kk)).contributor;
                    end;
                otherwise
                    disp('Unsupported subscripted reference');
                    disp_usage;
                    return;
            end;
            if length(ind)==1
                out = out{1};
            end;
        end;
    case 'origin'
        if isempty(myEvent.origin)
            beep;
            disp('Warning: Event informations are not retrieved yet!')
        end;
        if length(s)==1
            out=myEvent.origin;
        elseif length(s) == 2
            if ~(strcmp(s(2).type,'()') | strcmp(s(2).type,'{}'))
                disp('Unsupported subscripted expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension subscripting not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            for kk=1:length(ind)
                out{kk}=myEvent.origin{ind};
            end;
        end;
    otherwise
        disp('Unsupported subscripted expression');
        disp_usage;
        return;
end;

function disp_usage
help MatEvent\subsref.m

return;
