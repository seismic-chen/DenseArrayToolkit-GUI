function [stack, besth, bestk, Eellipse] = hkstack(SEIS,T0,DT,P,H,KAPPA,VP,W,ISMTH,ISADD)
%
% Stack the amplitudes of a set of receiver functions at predicted
% times of a Ps conversion and the first two multiples for different
% crustal thickness and Vp/Vs ratios, following the method of Zhu and
% Kanamori (2000), JGR.  The only difference to their code is that
% their is no option for smoothing here.  You can apply that to your
% RFn data beforehand using a filter of your choice
%
% IN:
% SEIS = array of receiver function amplitudes [nt x nrf]
% T0 = start time for all receiver functions (s)
% DT = sample interval for all receiver functions (s)
% P = slowness for each RFn [nrf]
% H = Array of equally spaced depth values to investigate (km)
% KAPPA = Array of arbitrarily spaced Vp/Vs ratios to investigate 
% VP = (optional) crustal Vp (km/s) to use, default = 6.3 km/s
% W = (optional) weights for Ps,  PpPs and PsPs+PpSs (default: 0.7/0.2/0.1)
% ISADD = if true, include Psx3, Pp(Psx3), Psx4, Psx5, Pp(Psx5), Psx6 
%         multiples (default = false) 
%
%        Copyright (c) 1996-2006, L. Zhu (lupei@eas.slu.edu)

%--- hkstack.m --- 
% 
% Filename: hkstack.m
% Description: See Above
% Author: Lupei Zhu, March, 1997 at Caltech
% Maintainer: Iain Bailey 
% Created: Tue Feb 15 15:54:13 2011 (-0800)
% Version: 1
% Last-Updated: Fri Jun 17 13:19:52 2011 (-0700)
%           By: Iain Bailey
%     Update #: 190
% Compatibility: Matlab R2009a
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%-- Commentary: 
% 
% Transcribed Lupei Zhu's code almost exactly. Took out loops so we
% can take advantage of matlab's matrix operations.  Also took out
% smoothing and parameter smtht since this can be applied to the
% receiver functions in advance. Also don't include variations of
% program for single h or single kappa values.
%
% Not certain what the am_cor variable comes from
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%-- Change Log:
% Dec. 8, 2015, Yunfeng Chen, Add output for varstack 
% Dec. 10, 2015, Yunfeng Chen, calculat the error ellipse
% Mar. 23, 2017, Yunfeng Chen, add t0p3s and t2p3s phases to stacking
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%-- Code:


% In the original program, smoothing can be applied during the
% stacking.  Here this is not included... it can be applied to the
% traces prior to input

% set default properties
if( nargin < 9 ),
  ISADD = false;
  if( nargin < 8 ),
    W = [ 0.7, 0.2, 0.1 ]; % Default weights for Ps PpPs PsPs/PpSs
    if( nargin < 7 ),
      VP = 6.3; % default Vp (km/s)
    end
  end
end

% get dimensions
nh = numel(H);  % number of depths to investigate
nk = numel(KAPPA); % number of kappa values
nrf = numel( P ); % number of rfns

% check the orientation of the seis array
if( size( SEIS,2 ) ~= nrf ),
  SEIS = SEIS.';
  if( size( SEIS,2 ) ~= nrf ),
    error(['One of the SEIS array dimensions [%i x %i] ',...
	   'should be the same as numel(P) [%i]'], ...
	  [size(SEIS,1), size(SEIS,2), nrf]);
    return
  end
end

nt = size( SEIS,1 ); % number of time samples

% check the orientation of P
if( size( P, 1 ) > size( P, 2) ), P = P.'; end

% set up the output
stack = zeros( nh, nk, nrf ); % results of the grid search
besth = H(1); % best values of h
bestk = KAPPA(1); % best Vp/Vs

% add smoothing to weights
W = W / (2*ISMTH);

% Velocity params
invVp2 = 1./VP^2;
invVs2 = invVp2 * (KAPPA'.^2);  % for all values of kappa 

% For all receiver functions...
p2 = P.*P;  % horizontal slowness squared
am_cor = 151.5478.*p2 + 3.2896.*P + 0.2618; % amp correction for incidence

% For all receiver functions and all values of kappa...
% columns for values of p, repeated rows for values of kappa to test
eta_p = repmat( sqrt( invVp2 - p2 ), nk,1);  % vertical slowness or sth like that
eta_s = sqrt( repmat(invVs2, 1,nrf) - repmat(p2, nk,1) ); % vertical slowness for s

% Get times for initial depths, all values of p and kappa in normalised time since
% beginning of trace
tPs = ( H(1).*( eta_s - eta_p ) - T0 )/DT + 1;
tPpPs = ( H(1).*(eta_s + eta_p ) - T0 )/DT + 1;
tPsPs = ( H(1).*2*eta_s - T0 )/DT + 1;

% get the change in normalised time by a step increase in depth (all p and kappa)
dh = H(2)-H(1);
dtPs = ( dh.*( eta_s - eta_p ) )/DT;
dtPpPs = ( dh.*(eta_s + eta_p ) )/DT;
dtPsPs = ( dh.*2*eta_s )/DT;

% loop through all values of crustal thickness
for hi = 1:nh,
  
  % loop through values of kappa
  for ki = 1:nk,
    
    % get the smoothed amplitude for each receiver function
    amp_tPs = amp( tPs(ki,:)-ISMTH, tPs(ki,:)+ISMTH, SEIS );
    amp_tPpPs = amp( tPpPs(ki,:)-ISMTH, tPpPs(ki,:)+ISMTH, SEIS );
    amp_tPsPs = amp( tPsPs(ki,:)-ISMTH, tPsPs(ki,:)+ISMTH, SEIS );

    % stack the different arrivals
    stack(hi,ki,:) =  W(1)*amp_tPs./am_cor ...
	+ W(2)*amp_tPpPs ...
	- W(3)*amp_tPsPs;    
    
    % Extra multiples 
    if (ISADD),
      amp_tPsPsPs = amp( tPs(ki,:)+tPsPs(ki,:)-ISMTH, ...
			 tPs(ki,:)+tPsPs(ki,:)+ISMTH, SEIS );
      amp_tPpPsPsPs = amp( tPpPs(ki,:)+tPsPs(ki,:)-ISMTH, ...
			   tPpPs(ki,:)+tPsPs(ki,:)+ISMTH, SEIS );
      amp_tPsPsPsPs = amp( 2.*tPsPs(ki,:)-ISMTH, 2.*tPsPs(ki,:)+ISMTH, SEIS );

      amp_tPsPsPsPsPs = amp( tPs(ki,:)+2.*tPsPs(ki,:)-ISMTH, ...
			     tPs(ki,:)+2.*tPsPs(ki,:)+ISMTH, SEIS );
     
      amp_tPpPsPsPsPsPs = amp( tPpPs(ki,:)+2.*tPsPs(ki,:)-ISMTH, ...
				tPpPs(ki,:)+2.*tPsPs(ki,:)+ISMTH, SEIS );
      amp_tPsPsPsPsPsPs = amp( 3.*tPsPs(ki,:)-ISMTH, 3.*tPsPs(ki,:)+ISMTH, SEIS );
      stack(hi,ki,:) = squeeze(stack(hi,ki,:))' ...
	  - W(1)*amp_tPsPsPs - W(2)*amp_tPpPsPsPs + W(3)*amp_tPsPsPsPs + ...
	  - W(1)*amp_tPsPsPsPsPs - W(2)*amp_tPpPsPsPsPsPs + W(3)*amp_tPsPsPsPsPsPs;
    end
  
  end

  % add on time for next depth step
  tPs = tPs + dtPs;
  tPpPs = tPpPs + dtPpPs;
  tPsPs = tPsPs + dtPsPs;

end

% get the standard deviation and mean
varstack = reshape( var( stack, 0, 3 ), nh,nk).';
stack = reshape( mean( stack, 3 ), nh,nk).' ;

% get the best combination based on max of the mean
[maxA, i] = max(max(stack));
[i,j] = find( stack == maxA ,1);
besth = H(j);
bestk = KAPPA(i);

% TODO estimate the confidence bounds for besth and bestk
% Yunfeng Chen, Dec. 10, 2015
% test data
% load('hk.xyz');
% load('std.xyz');
% H = unique(hk(:,1));
% KAPPA = unique(hk(:,2));
% V = hk(:,3);
% varV = std(:,3);
dh = H(2)-H(1);
dk = KAPPA(2)-KAPPA(1); 
nh = numel(H);
nk = numel(KAPPA);
% stack = [];
% for j = 1:nk
%     for i = 1:nh
%         k = i+(j-1)*nh;
%         stack(i,j) = V(k);
%         varstack(i,j) = varV(k);
%     end
% end
% zhu's code stores the results as negative number, depth is in first
% dimension
f = -stack';
varf = varstack';
[f0,i] = min(f(:));
% check to see if f0 is located at the boundary
[ix,iy] = ind2sub(size(f),i);
if ix <= 2 || ix >= nh-1 || iy <= 2 || iy >= nk-1
    % unable to calculate uncertainty
      Eellipse.hstd = 999;
      Eellipse.kstd = 999;
      Eellipse.X = [];
      Eellipse.Y = [];
else
    dhh = f(i+1)+f(i-1)-2*f0;
    dkk = f(i+nh)+f(i-nh)-2*f0;
    dhk = 0.25*(f(i+nh+1)+f(i-nh-1)-f(i+nh-1)-f(i-nh+1));
    eps = 4*dhh*dkk-dhk*dhk;
    [ih,ik] = find( f == f0 ,1);
    h = ih-(2*dkk*(f(i+1)-f(i-1))-dhk*(f(i+nh)-f(i-nh)))/eps;
    k = ik-(2*dhh*(f(i+nh)-f(i-nh))-dhk*(f(i+1)-f(i-1)))/eps;
    f0 = f0-0.5*(h-ih)*(h-ih)*dhh-0.5*(k-ik)*(k-ik)*dkk-(h-ih)*(k-ik)*dhk;
    vh = dhh/dh^2;
    vk = dkk/dk^2;
    vhk = -dhk/(dh*dk);
    dev = varf(ih,ik);
    minH = H(1);
    maxH = H(end);
    mink = KAPPA(1);
    maxk = KAPPA(end);
    % ratioH = 4/(maxH-minH);
    % ratiok = 2/(maxk-mink);
    % cxx = vh/(ratioH*ratioH);
    % cyy = vk/(ratiok*ratiok);
    % cxy = vhk/(ratioH*ratiok);
    cxx = vh;
    cyy = vk;
    cxy = vhk;
    
    s1 = 0.5*(cxx-cyy);
    s2 = sqrt(s1*s1+cxy*cxy);
    theta = 0.5*atan2(cxy,s1)*180/pi;
    s1 = 0.5*(cxx+cyy)+s2;
    s2 = 0.5*(cxx+cyy)-s2;
    % draw the ellipse
    a = sqrt(dev/s1);
    b = sqrt(dev/s2);
    a0 = minH+dh*(ih-1);
    b0 = mink+dk*(ik-1);
    angle = theta;
    steps = 36;
    beta = -angle*(pi/180);
    sinbeta = sin(beta);
    cosbeta = cos(beta);
    alpha = linspace(0, 360, steps)' .* (pi / 180);
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);
    
    X = a0 + (a * cosalpha * cosbeta - b * sinalpha * sinbeta);
    Y = b0 + (a * cosalpha * sinbeta + b * sinalpha * cosbeta);
    
    Eellipse.X = X;
    Eellipse.Y = Y;
    Eellipse.hstd = (max(X)-min(X))/2;
    Eellipse.kstd = (max(Y)-min(Y))/2;
end
% % plot HK
% imagesc(H,KAPPA,stack');hold on;



% plot(X,Y, 'r.-');hold on;
return 

%------------------------------------------------------------
function a = amp( t1, t2, rfseis )
% integrate data[n] between t1 and t2 (normalized time by dt)
%
% t1 = start (normalised) times for each trace
% t2 = end (normalised) times for each trace
% rfseis = matrix of column seismograms

nt = size(rfseis,1);
nseis = size(rfseis,2);
am = zeros(1,nseis);
t1 = max(1,t1) + nt.*((1:nseis)-1); % convert to positions in the matrix
t2 = min(nt,t2) + nt.*((1:nseis)-1);

it1 = floor(t1); % get the next lowest integer to t1
i = it1 + 1;
dd = i-t1; % 1 - remainder of t1
am = dd.*( dd.*rfseis(it1) + (2.0-dd).*rfseis(i) );

it2 = ceil(t2); % integer higher than t2
for k=1:nseis,
  am(k) = am(k) + sum(rfseis(i(k):(it2(k)-1))) + sum(rfseis((i(k)+1):it2(k)));
end
dd = it2-t2;
am = am - dd.*(dd.*rfseis(it2-1)+(2.-dd).*rfseis(it2));

a = 0.5*am;

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         Copyright (c) 1996-2006, L. Zhu (lupei@eas.slu.edu)
% 
% Permission to use, copy, modify, and distribute this package and supporting
% documentation for any purpose without fee is hereby granted, provided
% that the above copyright notice appear in all copies, that both that
% copyright notice and this permission notice appear in supporting
% documentation.
% 
% In case that you modify any codes in this package, you have to
% change the name of the modified code. You are welcome to send me a note
% about your modification and any suggestion.
% 
% In case that you redistribute this package to others, please send me
% the contacting info (email addr. preferred) so that future updates
% of this package can be sent to them.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% hkstack.m ends here
