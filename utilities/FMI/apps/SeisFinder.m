function [event, seiswav, seisatt,channelinfo] = SeisFinder(eventServer, seisServer)

% SeisFinder: a function to get seismograms and their attribution
% as well as their channel info for one selected event. 
%
% Outputs:
%       event: a structure array containing the tracking event 
%              information
%       seiswav: a cell array containing retrieved seismograms
%       seisatt: a structure array saved attributions of retrieved
%                seismograms
%       channelinfo: a structure array kept channel information
%
% Attention: 
% For multiple retrieves, the function return just reflects the
% last part. All retrieved seismograms and channel information
% are saved into a set of files named by user. If all these files
% can be merge into one file without any problem by using "SAVE"
% button, the function will return all retrieved information
%
% Also check: MatEvent, MatSeismogram, get

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Sep, 2003

%close all;
%clear all;

if exist('eventServer')~=1
    eventServer = 'iris';
end;
if exist('seisServer')~=1
    seisServer = 'pond';
end;

if exist('myEvent')~=1
    myEvent=MatEvent(eventServer);
end;
mySeis = MatSeismogram(seisServer);
[event, seiswav, seisatt,channelinfo]=seis_Find(myEvent,mySeis);
clear myEvent mySeis;
