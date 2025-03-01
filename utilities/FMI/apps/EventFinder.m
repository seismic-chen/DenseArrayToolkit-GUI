function [events,param, selection] = EventFinder(eventSrever);

% EventFinder: a program to search all available events on IRIS
%   given a set of parameters, such as earthquake location, magnitude,
%   depth, etc. It returns three outputs. The first one is param,
%   which saves parameters one input for searching earthquakes.
%   The second one is a MatEvent object which save all found events
%   The third one is selection saving index numbers of events one
%   selected.
%
%   also check: MatEvent

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Aug, 2003

%clear all
%close all

if exist('eventServer')~=1
    eventServer = 'iris';
end;

if exist('myEvent')~=1
    myEvent=MatEvent(eventServer);
end;
[param, events, selection]=event_Find(myEvent);
clear myEvent;
