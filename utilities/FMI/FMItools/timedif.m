function td = timedif(t1, t2)

% TIMEDIFF calculate time difference in seconds between two time instance
%
% timediff(t1,t2) - calculate time difference in seconds between t1 and t2
% t1 and t2 may be FMI time string or vector format. 
%
% If t2 is later than t1, then it returns a positive value. Otherwise it
% returns a negative value
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also FMITime local2UTC string2time time2string timecmp timeshift
%

if isa(t1,'char')
    t1=string2time(t1);
end;
if isa(t2,'char')
    t2=string2time(t2);
end;

td = etime(t2, t1);
