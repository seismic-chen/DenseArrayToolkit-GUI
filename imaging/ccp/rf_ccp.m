function [cp, RayMatrix0, MidPoints0] = rf_ccp( p, backaz, dz, zmax, z, vp, vs, lat, lon, model_type)
% calculate the PmS conversion point of RF for a given velocity model
% Input: 
%       p: ray parameter in s/km
%       backaz: back-azimuth
%       dz: depth increment
%       zmax: maximum depth
%       vp: P-wave velocity
%       vs: S-wave velocity
%       lat: station latitude
%       lon: station longitude
%       model_type: flat or spherical earth, default is flat 
% Output:
%       cp: struct contains location of conversion point
%       RayMatrix0: a nz*nx*7 ray piercing point matrix where the columns are
%       RRF Amplitudes, depth (these two varaibles will be saved to matrix later),
%       Ps - lat,lon, and Pp - lat,lon and index of receiver function (from 1 to nx)
%       MidPoints0: a (nz-1)*nx*7 matrix where the columns are Ps - lon,lat,raylength
%       and Pp - lon,lat,raylength, index of receiver function (from 1 to nx).
% Feb. 9, 2023, Yunfeng Chen, Global Seismology Group, Zhejiang University
if nargin == 9
    model_type = 'flat';
end

EPS = 1e-6;
R = 6371; % radius of earth
nx = numel( p );
% force rayp to be row vector
if( nx > 1 )
  p=p(:)';
end

% get the depths
zpos = (0.0:dz:zmax)';
nz = numel(zpos);
% deal with discontinuities in the vel model
idisc = find( z(1:end-1) == z(2:end) );
z(idisc) = z(idisc) - EPS;

r = repmat(R - zpos(1:end-1),1,nx);

% interpolate the vel model between each layer
vp = interp1( z, vp, zpos(1:end-1)+0.5*dz, 'linear','extrap');
vs = interp1( z, vs, zpos(1:end-1)+0.5*dz, 'linear','extrap');

% repeat matrices
p = repmat( p, nz-1, 1 );
% convert rayp to s/rad
if strcmp(model_type,'spherical')
    p = p*R;
end

backaz = repmat( backaz, nz-1, 1 );
vp = repmat( vp, 1, nx );
vs = repmat( vs, 1, nx );
% get associated vertical slowness
switch model_type
    case 'flat'
        % get horizontal position
        qa = sqrt(1./vp.^2 - p.^2);
        qb = sqrt(1./vs.^2 - p.^2);
        dha = p.*dz./qa;
        dhb = p.*dz./qb;
        dla=dz./(qa.*vp);
        dlb=dz./(qb.*vs);
    case 'spherical'
        qa = sqrt((r./vp).^2 - p.^2);
        qb = sqrt((r./vs).^2 - p.^2);
        dha = p.*dz./(r.*qa);
        dhb = p.*dz./(r.*qb);
        dla=r.*dz./(qa.*vp);
        dlb=r.*dz./(qb.*vs);
        % convert rad to degree
        dha = dha .* (360/(2*pi));
        dhb = dhb .* (360/(2*pi));
        % convert to km
        dha = dha .* (2*pi*R)/360;
        dhb = dhb .* (2*pi*R)/360;
end

% search for imaginary values in the distance calculations
for n = 1:size(dha,2)
    StopIndex = find(imag(dha(:,n)),1);
    if ~isempty(StopIndex)
        dha(StopIndex:end,n) = NaN * ones(length(dha(:,n))-StopIndex+1,1);
        dhb(StopIndex:end,n) = NaN * ones(length(dhb(:,n))-StopIndex+1,1);
    end
end

hposa = cumsum( dha, 1 );   % horizontal position of P phase piercing point
eposa = [zeros(1,nx) ; hposa.*sind(backaz)]; % easting of P phase piercing point
nposa = [zeros(1,nx) ; hposa.*cosd(backaz)]; % northing of P phase piercing point
la = dla;  % length of raypath segment of P phase

hposb = cumsum( dhb, 1 );   % horizontal position of Ps phase piercing point
% hposb = cumsum((dz./r) ./ sqrt(1./(p^2.*(r./vs).^-2))-1);
eposb = [zeros(1,nx) ; hposb.*sind(backaz)]; % easting of Ps phase piercing point
nposb = [zeros(1,nx) ; hposb.*cosd(backaz)]; % northing of Ps phase piercing point
lb = dlb;  % length of raypath segment of P phase

% calculate the lat lon of the conversion point
arclena = km2deg(hposa);
arclenb = km2deg(hposb);
az = backaz;

lat = repmat(lat, nz-1, 1);   % latitude of station
lon = repmat(lon, nz-1, 1);   % longitude of station
[lata,lona] = reckon(lat,lon,arclena,az);
[latb,lonb] = reckon(lat,lon,arclenb,az);
% output conversion point information as a structure array
cp=[];
for n=1:nx
    tmp.lata = lata(:,n);
    tmp.lona = lona(:,n);
    tmp.latb = latb(:,n);
    tmp.lonb = lonb(:,n);
    tmp.zpos = zpos(2:end);
    tmp.hposa = hposa(:,n);
    tmp.nposa = nposa(:,n);
    tmp.eposa = eposa(:,n);
    tmp.la=la(:,n);
    tmp.hposb = hposb(:,n);
    tmp.nposb = nposb(:,n);
    tmp.eposb = eposb(:,n);
    tmp.lb=lb(:,n);
    cp = [cp, tmp];
end
%% Get the RayMatrix and MidPoints matrices
nz=length(zpos);
RayMatrix0=NaN*ones(length(zpos),nx,7);
RayMatrix0(2:end,:,3)=[cp.latb];
RayMatrix0(2:end,:,4)=[cp.lonb];
RayMatrix0(2:end,:,5)=[cp.lata];
RayMatrix0(2:end,:,6)=[cp.lona];
% now deal with the first point, which is the station location
RayMatrix0(1,:,3)=lat(1,:);
RayMatrix0(1,:,4)=lon(1,:);
RayMatrix0(1,:,5)=lat(1,:);
RayMatrix0(1,:,6)=lon(1,:);
RayMatrix0(1:end,:,7)=ones(nz,1)*[1:nx];
% compute the mid points, which is roughly the average of each ray segment
MidPoints0=0.*ones(length(zpos)-1,nx,7);
MidPoints0(:,:,3)=[cp.lb];  % length of raypath segment of Ps
MidPoints0(:,:,6)=[cp.la];  % length of raypath segment of P
for n=1:nz-1
   MidPoints0(n,:,1)=(RayMatrix0(n,:,3)+RayMatrix0(n+1,:,3))/2;  % Ps midpoint lat
   MidPoints0(n,:,2)=(RayMatrix0(n,:,4)+RayMatrix0(n+1,:,4))/2;  % Ps midpoint lon
   MidPoints0(n,:,4)=(RayMatrix0(n,:,5)+RayMatrix0(n+1,:,5))/2;  % P midpoint lon
   MidPoints0(n,:,5)=(RayMatrix0(n,:,6)+RayMatrix0(n+1,:,6))/2;  % P midpoint lon
   MidPoints0(n,:,7)=RayMatrix0(n,:,7);  % Index of RF
end