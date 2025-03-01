
% FMI (Fissure-to-Matlab Interface) time format
%
% FMI supports two types of time representation, string and vector. All
% time used in FMI is UTC. If your time is local time, use local2UTC 
% function to convert. 
%
% The time string in FMI is a simplified version of ISO 8601 standard.
% Following time formats are accepted. Actually these time formats are same
% as the formats used in original FISSURE, expect the support to timezone
% is removed for simplification.
%
%   yyyy-mm-ddThh:mm:ss.dddZ 
%   yyyy-mm-ddThhmmss.dddZ 
%   yyyymmddThh:mm:ss.dddZ 
%   yyyymmddThhmmss.dddZ 
%   yyyyjjjThh:mm:ss.dddZ
%   yyyyjjjThhmmss.dddZ
%
%   yyyy - year
%   mm - month
%   jjj - day of the year
%   dd - day
%   T - seperator
%   hh - hour
%   mm - minutes
%   ss - seconds
%   ddd - milliseconds.
%   Z - indicates the time is UTC (Universal Time)
%
% The time vector is a six-element row vector [year, month, day, hour, minute, second].
% Except the second, all other fields must be integer. 
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   Sept, 2003
%
% See also local2UTC time2string string2time timecmp timeshift timedif
%

help FMITime
