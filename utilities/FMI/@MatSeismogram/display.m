function display(mySeis)

% DISPLAY Class MatNetwork display method
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

disp(' ')
disp('Class MatSeismogram');
disp(sprintf('  %d seismograms available',length(mySeis.availSeisFilt)));
disp(sprintf('  %d seismograms selected',length(mySeis.selectedSeis)));
disp(sprintf('  %d seismograms retrieved',mySeis.nSeismogram));
disp(' ');
