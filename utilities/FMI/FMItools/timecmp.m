function r=timecmp(t1,t2)

% TIMECMP compare two times.
% 
% r=timecmp(t1,t2) returns 1 if t1 is later than t2, -1 if t1 is earlier 
% than t2, 0 if same. t1 and t2 could be either string or vector
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also FMITime local2UTC time2string string2time timeshift timedif
%

if isa(t1,'char')
    t1=string2time(t1);
end;
if isa(t2,'char')
    t2=string2time(t2);
end;

tn1=t1(1)*1e8 + t1(2)*1e6 + t1(3)*1e4 + t1(4)*1e2 + t1(5);
tn2=t2(1)*1e8 + t2(2)*1e6 + t2(3)*1e4 + t2(4)*1e2 + t2(5);

if tn1>tn2
    r=1;
end;
if tn1<tn2
    r=-1;
end;
if tn1==tn2
    sec1 = t1(6);
    sec2 = t2(6);
    if sec1==sec2
        r = 0;
    elseif sec1>sec2
        r = 1;
    else
        r = -1;
    end;
end;
