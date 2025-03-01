function out = get(myEvent,name)

% get(myEvent,'PropertyName') class MatEvent subscripting expression interpreter
%
% Following property names are acceptable:
%   "name"              - Names of events
%   "nEvent"            - Number of retrievd events
%   "origin"            - Origins of events
%   "preferredOrigin"   - Preferred origins of events
%  
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

out = [];

switch lower(name)
    case 'nevent'
        out=myEvent.nEvent;
    case 'name'
        out=myEvent.name;
    case 'origin'
        out=myEvent.origin;
    case 'preferredorigin'
        out=myEvent.preferredOrigin;
    otherwise
        disp('Unrecognized property name');
        help MatEvent\get.m;
        return;
end;
