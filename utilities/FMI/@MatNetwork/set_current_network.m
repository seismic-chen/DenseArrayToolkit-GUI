function myNetwork=set_current_network(myNetwork,n)

% SETCURRENTNETWORK  set field currentNetwork 
%
% myNetwork=set_current_network(myNetwork,index)
% 
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

myNetwork.currentNetwork=n(1);
myNetwork.networkStation=struct([]);
myNetwork.findStationDone=0;
