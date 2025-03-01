function contributors = known_contributor(myEvent)

% KNOWN_CONTRIBUTOR find known catalogs
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

disp('  Retrieving information for DMC, Please wait ... ');
contributors = cell(myEvent.eventFnd.known_contributors);
