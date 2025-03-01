function display(myNetwork)

% DISPLAY Class MatNetwork display method
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

disp(' ')
disp('Class NetworkDC');
disp(sprintf('   %d networks found',myNetwork.networkNum));
disp(sprintf('   %d stations found in current chosen network %s', ... 
        myNetwork.networkStationNum, myNetwork.networkInfo(myNetwork.currentNetwork).code) );
disp(sprintf('   %d stations found in station pool',myNetwork.stationPoolNum));
disp(' ');
