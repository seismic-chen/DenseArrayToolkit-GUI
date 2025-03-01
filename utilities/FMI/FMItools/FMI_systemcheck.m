function status_out = FMI_systemcheck;

% FMI system requirement and setting check
% 
% If everything is ok, FMI_systemcheck will return 1; otherwise return 0.
% 
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

% check FMI path
p = version;
if strcmp(p(1:3),'6.5') == 0;
    status = 0;
    disp('The matlab version must be 6.5 (R13) or above');
    return;
end;

status = 1;
error = 0;
p=which('MatEvent');
if isempty(p)
    error = 1;
    disp('Cannot locate MatEvent class.');
end;
p=which('MatNetwork');
if isempty(p)
    error = 1;
    disp('Cannot locate MatNetwork class.');
end;
p=which('MatSeismogram');
if isempty(p)
    error = 1;
    disp('Cannot locate MatSeismogram class.');
end;
p=which('FMITime');
if isempty(p)
    error = 1;
    disp('Cannot locate FMItools directory.');
end;
p=which('m_coast');
if isempty(p)
    error = 1;
    disp('Cannot locate M_Map directory.');
end;
p=which('taupTime');
if isempty(p)
    error = 1;
    disp('Cannot locate MatTaup directory.');
end;

if error
    disp('Please check if FMI paths have been properly set');
    status = 0;
else
    %disp('FMI paths OK!');
end;

% check Java dependencies
error = 0;
s=methods('edu.iris.Fissures.IfEvent.EventFinder');
if isempty(s)
    error = 1;
    disp('Cannot locat fissuresIDL');
end;
s=methods('edu.iris.Fissures.event.EventAttrImpl');
if isempty(s)
    error = 1;
    disp('Cannot locate fissuresImpl');
end;
s=methods('edu.sc.seis.fissuresUtil.namingService.FissuresNamingService');
if isempty(s)
    error = 1;
    disp('Cannot locate fissureUtil');
end;
s=methods('junit.framework.Assert');
if isempty(s)
    error = 1;
    disp('Cannot locate junit');
end;
s=methods('junitx.framework.Assert');
if isempty(s)
    error = 1;
    disp('Cannot locate junit-addon');
end;
s=methods('org.omg.CORBA_2_3.ORB');
if isempty(s)
    error = 1;
    disp('Cannot locate OB');
end;
s=methods('com.ooc.CosNaming.OBNamingContext');
if isempty(s)
    error = 1;
    disp('Cannot locate OBNaming');
end;
s=methods('edu.iris.dmc.seedcodec.Codec');
if isempty(s)
    error = 1;
    disp('Cannot locate SeedCodec');
end;
s=methods('edu.sc.seis.TauP.MatTauP_Time');
if isempty(s)
    error = 1;
    disp('Cannot locate matTaup');
end;
s=methods('edu.sc.seis.TauP.TauP_Time');
if isempty(s)
    error = 1;
    disp('Cannot locate TauP');
end;

if error
    disp('Please check if FMI class paths have been properly set');
    status = 0;
else
    %disp('FMI Java class paths OK!');
end;

if nargout == 0 
    if status
        disp('FMI is installed successfully!'); disp(' ');
    else
        disp('FMI is not installed properly!');
    end;
else
    status_out = status;
end;
