function [timeout, seisout, depth0] = rf_migrate( time, seis, p, dz, zmax, z, vp, vs, TimeCorrections)
% conduct the time-to-depth conversion of receiver functions
% Input: 
%       time: cell array contains time axis for all traces
%       seis: cell array contains waveforms for all traces
%       p: vector contains all ray parameters in s/km
%       dz: depth increment
%       zmax: maximum depth
%       vp: P-wave velocity
%       vs: S-wave velocity
%       TimeCorrections: amount of time corrections
% Output:
%       timeout: cell array contains the output time axis
%       seisout: cell array contains the migrated waveforms
%       depth0:  depth axis
% Feb. 9, 2023, Yunfeng Chen, Global Seismology Group, Zhejiang University
EPS = 1e-6;

nseis = numel( p );

% get the depths
zout = (0.0:dz:zmax);
% check if TimeCorrections term exist
if nargin == 8 
    TimeCorrections = zeros(length(zout),length(seis)); 
end

% deal with discontinuities in the vel model
idisc = find( z(1:end-1) == z(2:end) );
z(idisc) = z(idisc) - EPS;

% interpolate the vel model in middle of each interval
vp = interp1( z, vp, zout+0.5*dz, 'linear','extrap');
vs = interp1( z, vs, zout+0.5*dz, 'linear','extrap');
r = 6371-zout;
% convert unit from s/km to s/rad;
p = skm2srad(p);
% calculate the tPs for each RF
for iseis = 1:nseis
    % get vertical slowness for each depth point at this rayp
    qa = sqrt((r./vp).^2 - p(iseis).^2);
    qb = sqrt((r./vs).^2 - p(iseis).^2);
    % get time difference
    dt = (dz./r).*(qb - qa); % time difference at each depth interval
    
    % get time at bottom of each layer
    tout = cumsum(dt) ;
    % add time correction to Tpds differntial time calucated using 1D
    % velcoity
    tout = tout + TimeCorrections(:,iseis)';
    
    % get time at top of each layer
    % tout = tout - dt;
    
    % get rid of times and depths beyond the end of the seismogram
    %idx = find( tout < max( time{iseis} ) );
    %tout = tout( idx );
    
    tout( imag(tout)~=0 ) = Inf;   % correct for evernecent waves
    
    t = time{iseis};
    x = seis{iseis};
    xout = interp1( t, x, tout , 'linear', 0);
    % deal with the negaitve time
    if size(xout,2) > size(xout,1)
        xout = xout';
    end
    seisout{iseis} = xout;
    seisout{iseis}(tout == Inf) = NaN;   % correct for evernecent waves
    timeout{iseis} = tout;
    depth0 = zout;
end
return
