% Jan, 15,2017, Yunfeng Chen, Global Seismology Group, University of
% Alberta
function [best_pair] = hk_stack_bootstrap(stack1,stack2,stack3)
% compute the HK stacking given individual results for three phases and
% returns the best H and K pari
% Input: 
% stack1: a N by NH*NK matrix contains the amplitude of Pms phase at the
% predicated arrival time
% stack2: a N by NH*NK matrix contains the amplitude of PpPms phase at the
% predicated phase time
% stack3: a N by NH*NK matrix contains the amplitude of PpSms phase at the
% predicated arrival time
% Output:
% best_pair: 1*2 matrix containts the best H and K estimate corresponds 
% to the maximum amplidute in HK domain

% Hall = 30:0.1:70;
% % make sure the kappa is the same with that in the main function
% kappa = 1.6:0.01:2.0;
% kappa = 1.6:0.002:2.0;

global Hall
global kappa

[N,M] = size(stack1);
% Pms
Cn1 = sum(stack1,1);
% PpPms
Cn2 = sum(stack2,1);
% PpSms
Cn3 = sum(stack3,1);
w1 = 1/3;
w2 = 1/3;
w3 = 1/3;
Cn = [w1.*Cn1+w2.*Cn2+w3.*Cn3];
[~,index] = max(Cn);
[indexh,indexk] = ind2sub([length(Hall),length(kappa)],index);
besth = Hall(indexh);
bestk = kappa(indexk);
best_pair = [besth,bestk];