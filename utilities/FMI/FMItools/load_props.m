function props = load_props(nameServer);

% load_props loads naming server properties. The choices for naming servers
% are the IRIS DMC and USC SCEPP. Other servers may exist. Port 6371 are used 
% by both USC and the DMC.
%
% Input argument:
%    nameServer:  'IRIS_DMC' or 'USC_SCEPP'. The IRIS DMC is the default
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   April, 2004
%

import java.io.*;
import java.lang.*;

if exist('nameServer') ~=1 || isempty(nameServer)
    nameServer = 'IRIS_DMC';
end;

props = System.getProperties;
props.setProperty('org.omg.CORBA.ORBClass', 'com.ooc.CORBA.ORB');
props.setProperty('org.omg.CORBA.ORBSingletonClass','com.ooc.CORBA.ORBSingleton');
switch nameServer
    case 'IRIS_DMC'
        props.setProperty('ooc.orb.service.NameService', ...
                  'corbaloc:iiop:dmc.iris.washington.edu:6371/NameService');
    case 'USC_SCEPP'
        props.setProperty('ooc.orb.service.NameService', ...
                  'corbaloc:iiop:scepp.seis.sc.edu:6371/NameService');
    otherwise
        error('Unrecognized naming server');      
end
props.setProperty('ooc.orb.oa.conc_model', 'thread_pool');
props.setProperty('ooc.orb.oa.thread_pool', '5'); 
