function viewComponents(zfname)

% Given the sac filename for the z component, find the other
%  components (assuming same directory and BHZ, BHE, BHN suffixes),
%  then plot all
%
format compact;

% % add programs used to your path
% addpath ~/seismology/programs/processRFmatlab/ioFunctions/
% addpath ~/seismology/programs/processRFmatlab/plotFunctions/

% [tz, z, atimes, labels] = getseis( zfname ); % get the zcompontent
  [tz, z, hdrz] = sac2mat( zfname );
% find the corresponding e component
m = regexp( zfname,'.BHZ', 'split');
% efname = [ m{1}{1:1},'.BHE'];
% [te, e] = getseis( efname ); % get the zcompontent
  efname = [ m{1},'.BHE',m{2}];
  [te, e, hdre] = sac2mat( efname );
% and n component
% nfname = [ m{1}{1},'.BHN'];
% [tn, n] = getseis( nfname ); % get the zcompontent
  nfname = [ m{1},'.BHN',m{2}];
  [tn, n, hdrn] = sac2mat( nfname );
 

tmin = min( [tz, te, tn] );
tmax = max( [tz, te, tn] );
amin = min( [z', e', n'] );
amax = max( [z', e', n'] );
  
% plot
% plot3seis( tz, z, te, e, tn, n, [tmin, tmax, amin, amax], ...
% 	   atimes, labels);
  
plot3seis( tz, z, te, e, tn, n, [tmin, tmax, amin, amax]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [t, a, times, labels] = getseis( fname )
% 
% % from a file name, load the file and get the amplitude
% %disp(fname{:})
% [t, a, head] = sac2mat( fname );
% 
% % arrival times
% head.times.atimes
% 
% % times = [ head.times.atimes[t0 , head.times.atimes[t1 , ...
% % 	  head.times.atimes[t2 , head.times.atimes[t3 , ...
% % 	  head.times.atimes[t4 , head.times.atimes[t5 , ...
% % 	  head.times.atimes[t6 , head.times.atimes[t7 , ...
% % 	  head.times.atimes[t8 , head.times.atimes[t9  ...
% % 	];
% % 
% % % labels for arrival times
% % labels = strvcat( head.times.atimes[kt0 , head.times.atimes[kt1 , ...
% % 		  head.times.atimes[kt2 , head.times.kt3 , ...
% % 		  head.times.kt4 , head.times.kt5 , ...
% % 		  head.times.kt6 , head.times.kt7 , ...
% % 		  head.times.kt8 , head.times.kt9  ...
% % 		  );
% 
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
