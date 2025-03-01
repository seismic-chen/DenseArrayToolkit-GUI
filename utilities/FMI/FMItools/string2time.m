function time=string2time(s)

% STRING2TIME convert time string to time vector
%
% time=string2time(s)
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also FMITime local2UTC time2string timecmp timeshift timedif
%

if ~isnumeric(s) & ~ischar(s)
    error('Input time format must be a character array');
end;
if isnumeric(s)
    time = s;
    return;
end;

s=strrep(s,'-','');
s=strrep(s,':','');

index = find(s == 'T');
if isempty(index) || (index~=8 && index~=9)
    error('Wrong time format');
end;

time = zeros(1,6);
if index == 9
    time = sscanf(s,'%4d%2d%2dT%2d%2d%f',6)';
else
    v = sscanf(s,'%4d%3dT%2d%2d%f',5);
    n = datenum(v(1),1,1);
    n = n+v(2)-1;
    d = datevec(n);
    
    time = [d(1) d(2) d(3) v(3) v(4) v(5)];
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
