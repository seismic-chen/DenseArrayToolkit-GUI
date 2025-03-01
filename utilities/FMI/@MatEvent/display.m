function display(myEvent)
% Class MatEvent display methods
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

disp('Object of the class EventDC');
fprintf('  %d events found.',myEvent.nEvent);
if myEvent.nEvent>0 & isempty(myEvent.name)
    disp('  Event information has not been retrieved');
end
disp(' ');
disp(' ');
