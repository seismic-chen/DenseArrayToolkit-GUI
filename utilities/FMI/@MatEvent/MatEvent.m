function myMatEvent = MatEvent(server)

% MatEvent class constructor
%
% myMatEvent = MatEvent(server) creates an object of MatEvent class
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

import edu.iris.Fissures.IfSeismogramDC.*;
import edu.iris.Fissures.IfNetwork.*;
import edu.iris.Fissures.IfEvent.*;
import edu.iris.Fissures.*;
import edu.iris.Fissures.model.*;
import edu.iris.Fissures.utility.*;
import edu.sc.seis.fissuresUtil.namingService.*;
import java.io.*;
import java.lang.*;

if exist('server')~=1 || isempty(server)
    server = 'iris';
end;

switch lower(server)
    case {'iris', 'iris_eventdc'}
        serverName = 'IRIS_EventDC';
        serverDNS = 'edu/iris/dmc';
    case {'berkeley', 'ncedc','ncedc_eventdc'}
        serverName = 'NCEDC_EventDC';
        serverDNS = 'edu/berkeley/geo';
    case {'berkeley-ncsn', 'ncsn','ncsn_eventdc'}
        serverName = 'NCSN_EventDC';
        serverDNS = 'edu/berkeley/geo';
    case {'usc', 'scepp', 'scepp_eventdc'}
        serverName = 'SCEPPEventDC';
        serverDNS = 'edu/sc/seis';
    case {'caltech','scedc','scedc_eventdc'}
        serverName = 'SCEDC_EventDC';
        serverDNS = 'edu/caltech/gps/k2';
    otherwise
        error('Unrcognized server');
end;

% load DMC naming server properties
props = load_props('IRIS_DMC');

try
    fissNaming=FissuresNamingService(props);
    %fissNaming=FissuresNamingServiceNoLogger(props);
    %disp('Naming service initialzed'); 
catch
    disp('Error: Cannot initialize Fissures naming service. The CORBA naming server');
    disp('may be down, or more likely there is a network problem for your computer');
    disp('or a firewall is blocking MATLAB from the Internet. Please check your');
    disp('network configuration and try again. If the problem still exists, report');
    disp('this problem to fmihelp@ess.washington.edu. You may also check status of');
    disp('the Fissures servers at http://seis.sc.edu/wily/Check');
    error('Cannot Fissures initialize naming service');
end;

try
    eventDC = fissNaming.getEventDC(serverDNS, serverName);
    disp('Got EventDC');
catch
    disp('Error: Cannot resolve EventDC. The event server may be down. Please');
    disp('report this problem to fmihelp@ess.washington.edu. You may also check');
    disp('status of the Fissures servers at http://seis.sc.edu/wily/Check');
    error('Cannot resolve EventDC');
end;
eventFnd = eventDC.a_finder;
%eventChanFnd = eventDC.a_channel_finder;

myMatEvent.eventFnd = eventFnd;
%myMatEvent.eventChanFnd = eventChanFnd;
myMatEvent.eventAccess=[];
myMatEvent.nEvent=0;
myMatEvent.name = {};
myMatEvent.origin = {};
myMatEvent.preferredOrigin = struct([]);
myMatEvent=class(myMatEvent,'MatEvent');

return;
