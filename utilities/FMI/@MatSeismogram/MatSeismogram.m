function myMatSeis=MatSeismogram(server)

% MATSEISMOGRAM class constructor
%
% myMatSeis = MatSeismogram(server) creates an object of MatSeismogram
% class. Server could be "Archieve", "Bud", or "Pond". The default is
% "Pond"
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
import edu.iris.Fissures.*;
import java.io.*;
import java.lang.*;
import edu.sc.seis.fissuresUtil.namingService.*;

if exist('server')~=1 || isempty(server)
    server = 'pond';
end;

switch lower(server)
    case {'pond', 'iris_pond', 'iris_ponddatacenter'}
        seisServerName = 'IRIS_PondDataCenter';
        netServerName = 'IRIS_NetworkDC';
        serverDNS = 'edu/iris/dmc';
    case {'bud', 'iris_bud', 'iris_buddatacenter'}
        seisServerName = 'IRIS_BudDataCenter';
        netServerName = 'IRIS_NetworkDC';
        serverDNS = 'edu/iris/dmc';
    case {'archive', 'iris_archive', 'iris_archivedatacenter'}
        seisServerName = 'IRIS_ArchiveDataCenter';
        netServerName = 'IRIS_NetworkDC';
        serverDNS = 'edu/iris/dmc';
    case {'caltech', 'scedc', 'scedc_datacenter'}
        seisServerName = 'SCEDC_DataCenter';
        netServerName = 'SCEDC_NetworkDC';
        serverDNS = 'edu/caltech/gps/k2';
    case {'berkeley', 'ncedc','ncedc_datacenter'}
        seisServerName = 'NCEDC_DataCenter';
        netServerName = 'NCEDC_NetworkDC';
        serverDNS = 'edu/berkeley/geo';
    case {'usc', 'scepp', 'scepp_datacenter'}
        seisServerName = 'SCEPPSeismogramDC';
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
    disp('the Fissures servers at http://piglet.seis.sc.edu:8080/wily/Check');
    error('Cannot initialize Fissures naming service');
end;

try
    netDc = fissNaming.getNetworkDC(serverDNS,netServerName');
    disp('Got NetworkDC');
catch
    disp('Error: Cannot resolve NetworkDC. The network server may be down. Please');
    disp('report this problem to fmihelp@ess.washington.edu. You may also check');
    disp('status of the Fissures servers at http://piglet.seis.sc.edu:8080/wily/Check');
    error('Cannot resolve NetworkDC');
end;
netExp = netDc.a_explorer;
%netFnd = netDc.a_finder;

try
    seisDc = fissNaming.getSeismogramDC(serverDNS,seisServerName);
    disp('Got SeismogramDC');
catch
    disp('Error: Cannot resolve SeismogramDC. The seismogram server may be down. Please');
    disp('report this problem to fmihelp@ess.washington.edu. You may also check');
    disp('status of the Fissures servers at http://piglet.seis.sc.edu:8080/wily/Check');
    error('Cannot resolve SeismogramDC');
end;

myMatSeis.seisDc = seisDc;
myMatSeis.netExp = netExp;
myMatSeis.server = upper(server);
myMatSeis.availSeisFilt = [];
% myMatSeis.availFrom = {};
myMatSeis.availChanLoc = [];

myMatSeis.selectedSeis = [];

myMatSeis.seismogram = [];
myMatSeis.seismogramAttr = struct([]);
myMatSeis.seisChanId = [];
myMatSeis.seisChanInfo = [];
myMatSeis.nSeismogram = 0;

iterator.enabled = 0;
iterator.nTotal = 0;
iterator.nPerRetrieve = 0;
iterator.nRetrieved = 0;
iterator.nRemained = 0;
iterator.nIter = 0;
myMatSeis.iterator = iterator;

myMatSeis=class(myMatSeis,'MatSeismogram');

return;
