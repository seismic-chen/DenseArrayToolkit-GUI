function [SEIS, KEY] = rotateSeisENZ2TLQ( SEIS , BAZ, INC )
%
%
% Rotate channels from E-N-Up to 
% Transverse-L(Along_Ray/P)-Q(Ray_Perpendicular/SV) coordinate system
%
% IN:
% SEIS = Seismogram array (NT x 3), 1 column for each component,
%        must have correct order [ T, R, Z ]  
% BAZ = back azimuth to eqk from station in degrees
% INC = Incident angle in degrees
%
% OUT:
% SEIS = rotated seismogram in new coord system
% KEY = names of each component
% Change Log:
% June 29, 2018, Y.C., the rotation matrix used by the original code was completely wrong, use correct one instead
% get components
E=SEIS(:,1);
N=SEIS(:,2);
Z=SEIS(:,3);

%L =  Z*cosd(INC) - sind(INC)*( N*cosd(BAZ) - E*sind(BAZ) );
%Q = -Z*sind(INC) - cosd(INC)*( N*cosd(BAZ) - E*sind(BAZ) );
%T = -N*sind(BAZ) + E*cosd(BAZ);

L =  cosd(INC)*Z - sind(INC)*sind(BAZ)*E - sind(INC)*cosd(BAZ)*N ;
Q =  sind(INC)*Z + cosd(INC)*sind(BAZ)*E - cosd(INC)*cosd(BAZ)*N ;
T = -cosd(BAZ)*E + sind(BAZ)*N;

SEIS = [ T, L, Q ];
KEY = strvcat('T','L','Q');
return;
