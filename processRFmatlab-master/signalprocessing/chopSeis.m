function [SEIS, T] = chopSeis( SEIS, T, T1, T2 )
% 
% [SEIS, T] = chopSeis( SEIS, T, T1, T2 )
% 
% Chop the seismograms so only T1 <= t <= T2 are included
%
% In: 
% SEIS = Seismogram array (NT x NC), 1 column for each component.  
% T = the times of the samples in  SEIS in s
% T1 = the new start time
% T2 = the new end time
%
% Out:
% SEIS = seismogram array after chop
% T = new time
%

% chopSeis.m --- 
%  
%  Filename: chopSeis.m
%  Author: Iain Bailey
%  Header Created: Thu Nov 17 15:04:23 2011 (-0800)
%  Version: 1
%  Last-Updated: Thu Nov 17 15:09:06 2011 (-0800)
%            By: Iain Bailey
%      Update #: 8
%  
%----------------------------------------------------------------------
%  
% Change Log:
%  Thu Nov 17 2011 : Corrected error reporting
%  
%  Apr. 22, 2015: Fill zero if signal is not long enough, Yunfeng Chen 
% Code:
% Jul 20, 2015: Check to see if the signal length is correct as the number
% of data points might not be exactly the same as (T2-T1)/dt dut to the
% precision issue when rounding the number.
% Mar. 5, 2016,: the end point related to T2 is determined by index of T1 +
% the length of signal, which garentee that the signal after cutting has
% the correct length (npts)
% Mar. 16, 2016: Force dt to keep three decimal place
% Apr. 16, 2016: Do not include the last point when choping the seismogram

if( T(1) > T1 ), 
  fprintf('Warning: Requested chop before %.2f, but time at trace start is %.2f\n',...
	  [T1, T(1)] );
end
if( T(end) < T2 ),
  fprintf('Warning: Requested chop after %.2f, but time at trace end is %.2f\n',...
	  [T2, T(end)] );
end
DT = round(T(2)-T(1),3);
% padding zeros if T(end) < T2 and T(1) > T1
if T(end) < T2
   nzero = ceil((T2-T(end))/(T(2)-T(1)));
   SEIS = padarray(SEIS,[nzero],'post');
   T = T(1):DT:T2;
end
if T(1) > T1
   nzero = ceil((T(1)-T1)/(T(2)-T(1)))
   SEIS = padarray(SEIS,[nzero],'pre');
   T = T1:DT:T(end);
end
% Yunfeng comment: only take the portion before T2 (Exclude the point that 
% is equal to T2)
% keepIdx = find( T >= T1 & T < T2 );
DT = T(2)-T(1);
npts = (T2-T1)/DT;
tmpIdx = int32(find( T >= T1 ));
startIdx = tmpIdx(1);
endIdx = startIdx + npts;
SEIS = SEIS( startIdx:endIdx-1, : );
T = T(startIdx:endIdx-1);
% % check to see if the signal length is correct
% if sum(keepIdx) ~= round((T2-T1)/(T(2)-T(1)))
%     num = sum(keepIdx)-round((T2-T1)/(T(2)-T(1)));
%     SEIS = SEIS(1:end-num,:);
%     T = T(1:end-num);
% end
return

%----------------------------------------------------------------------
%  
%  This program is free software; you can redistribute it and/or
%  modify it under the terms of the GNU General Public License as
%  published by the Free Software Foundation; either version 3, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; see the file COPYING.  If not, write to
%  the Free Software Foundation, Inc., 51 Franklin Street, Fifth
%  Floor, Boston, MA 02110-1301, USA.
%  
%----------------------------------------------------------------------
% chopSeis.m ends here
