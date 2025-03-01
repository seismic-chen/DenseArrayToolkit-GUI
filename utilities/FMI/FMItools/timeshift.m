function time_out=timeshift(time,shift)

% TIMESHIFT calculate time with a shift
%
% time_out=timeShift(time,shift)
% 
% Input Arguments:
%    time:  current time. It may be FMI time string or time vector;
%   shift:  shift is a number in seconds. It may be positive or negtive.
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also FMITime local2UTC time2string string2time timecmp timedif
%

if ischar(time)
    time=string2time(time);
end;

% check if date and time is valid.
if time(1) > 2100
    error('Wrong YEAR value');
end;
if time(2) > 12
    error('Wrong MONTH value');
end;
if time(3) > eomday(time(1),time(2))
    error('Wrong DAY value');
end;
if time(4) >= 24
    error('Wrong HOUR value');
end;
if time(5) >= 60 
    error('Wrong MINUTE value');
end;
if time(6) >= 60
    error('Wrong SECOND value');
end;

n = datenum(time(1:3));

limit = [0 0 0 24 60 60];
for ii = 6:-1:4 
    temp = time(ii)+shift;
    time(ii) = mod(temp,limit(ii));
    shift = (temp-time(ii))/limit(ii);
    if shift == 0
        break;
    end;
end;
if shift~=0
    n = n+shift;
end;

time_out = datevec(n);
time_out(4:6) = time(4:6);


