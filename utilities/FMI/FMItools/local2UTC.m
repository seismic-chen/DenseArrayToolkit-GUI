function UTC = local2UTC(localTime, td);

% LOCAL2UTC convert local time to UTC
%
% local2UTC(localTime, timeDiff) converts localTime, which may be string or
%   FMI time structure, to UTC by given time difference between local time
%   and UTC. timeDiff must be in the format of 'hh:mm' or 'hh', where hh is
%   difference in hours and mm is difference in minutes. If necessary, a
%   minus sign '-' should be applied.
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also FMITime time2string string2time timecmp timeshift timedif
%

if ischar(localTime)
    localTime=string2time(localTime);
end;

% check if date and time is valid.
if localTime(1) > 2100
    error('Wrong YEAR value');
end;
if localTime(2) > 12
    error('Wrong MONTH value');
end;
if localTime(3) > eomday(localTime(1),localTime(2))
    error('Wrong DAY value');
end;
if localTime(4) >= 24
    error('Wrong HOUR value');
end;
if localTime(5) >= 60 
    error('Wrong MINUTE value');
end;
if localTime(6) >= 60
    error('Wrong SECOND value');
end;

hh = sscanf(td,'%d:%d',2);
if isempty(hh)
    error('Wrong time difference format');
end;
if length(hh) == 1
    shift = hh*3600;
else
    if hh>=0
        shift = hh(1)*3600 + hh(2)*60;
    else
        shift = hh(1)*3600 - hh(2)*60;
    end;
end;

UTC = timeshift(localTime, shift);

