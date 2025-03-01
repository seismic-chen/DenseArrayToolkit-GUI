function SEIS = removeSeisTrend( SEIS )
% SEIS = removeSeisTrend( SEIS )
%
% Remove the trend in the seismogram
% 
% IN
% SEIS = Seismogram array (NT x NC), 1 column for each component. 
%
% Out:
% SEIS = seismograms after removing the trend
%
% Modified from IWB's codes by Yunfeng Chen, Feb. 2016

SEIS = detrend(SEIS);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
