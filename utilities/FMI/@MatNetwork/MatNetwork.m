function myNetworkDC = MatNetwork(server)

% MATNETWORK class constructor
%
% myNetworkDC = MatNetwork(server) creates an object of EventDC class
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

import edu.iris.Fissures.IfNetwork.*;
import edu.iris.Fissures.network.*;
import edu.iris.Fissures.*;
import edu.iris.Fissures.model.*;
import edu.iris.Fissures.utility.*;
import edu.sc.seis.fissuresUtil.namingService.*;
import java.io.*;
import java.lang.*;
import org.omg.CORBA.*;
import org.omg.CORBA.portable.*;
import org.omg.CosNaming.*;
import org.omg.CosNaming.NamingContextPackage.*;

if exist('server')~=1 || isempty(server)
    server = 'iris';
end;

switch lower(server)
    case {'iris', 'iris_networkdc'}
        netServerName = 'IRIS_NetworkDC';
        serverDNS = 'edu/iris/dmc';
    case {'caltech', 'scedc', 'scedc_networkdc'}
        netServerName = 'SCEDC_NetworkDC';
        serverDNS = 'edu/caltech/gps/k2';
    case {'berkeley', 'ncedc','ncedc_networkdc'}
        netServerName = 'NCEDC_NetworkDC';
        serverDNS = 'edu/berkeley/geo';
    case {'usc', 'scepp', 'scepp_networkdc'}
        netServerName = 'SCEPPNetworkDC';
        serverDNS = 'edu/sc/seis';
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
    error('Cannot initialize Fissures naming service');
end;

try
    netDc = fissNaming.getNetworkDC(serverDNS,netServerName);
    disp('Got NetworkDC');
catch
    disp('Error: Cannot resolve NetworkDC. The network server may be down. Please');
    disp('report this problem to fmihelp@ess.washington.edu. You may also check');
    disp('status of the Fissures servers at http://seis.sc.edu/wily/Check');
    error('Cannot resolve NetworkDC');
end;
netExp = netDc.a_explorer;
netFnd = netDc.a_finder;

myNetworkDC.networkFnd=netFnd;
myNetworkDC.networkExp=netExp;
myNetworkDC.networkInfo=struct([]);
myNetworkDC.networkAccess=[];
myNetworkDC.networkNum=0;
myNetworkDC.chosenNetwork=[];

myNetworkDC.networkStation=struct([]);
myNetworkDC.networkStationNum=0;
myNetworkDC.stationAddedToPool=0;
myNetworkDC.currentNetwork=[];
myNetworkDC.stationPool=struct([]);
myNetworkDC.stationPoolNum=0;
myNetworkDC.findStationDone=0;  % used if retrieve channels by station
myNetworkDC=class(myNetworkDC,'MatNetwork');

return;