function catalogs = known_catalog(myEvent)

% KNOWN_CATALOG find known catalogs
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

disp('  Retrieving information for DMC, Please wait ... ');
catalogs = cell(myEvent.eventFnd.known_catalogs);
