function s=time2string(time)

% TIME2STRING convert time structure to string
%
% s=time2string(time)
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also FMITime local2UTC string2time timecmp timeshift timedif
%

if ~isnumeric(time) & ~ischar(time)
    error('Input time format must be a 6-element numeric array');
end;
if ischar(time)
    s = time;
    return;
end;

% check if date and time is valid.
if time(1) > 3500
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

s=sprintf('%4d-%02d-%02dT%02d:%02d:%5.3fZ',time);
