% Calculate the average crustal velocity model using crust 1.0 model
% Nov 23, 2016, Yunfeng Chen, Global Seismology Group, University of
% Alberta.
% Use crust 1.0 model for HK calculation, note that the ray tracing code
% uses the top boundary layer conversion, means the velocity for each layer
% is define at top boundary. This is the same as crust 1.0 model
% read in crustal model
% Nov. 25, 2016, Bug fixed regarding syntehtics, the problem is the input
% format of the model, the format requirement is very stringent!
% Dec. 14, 2016, return the number of sedimentary layers, which is
% required by the HK program
% Dec. 19, 2016, using vp and vs if zero is not always good criterion to 
% count the number of sediments, I should also consider the layer thickness
% Jan. 25, 2017, add flag to control if include surface topography 
function varargout = obtain_crust1(varargin)
lat = varargin{1};
lon = varargin{2};
if nargin < 3
    H = 0;
else
    H = varargin{3};
end

if nargin < 4
    if_topo = 0;
else
    if_topo = varargin{4};
end

% lat = 54.5;
% lon = -115.5;
nla=180;
nlo=360;
nd=9;
mdldir='./velocity_model/crust1.0/';
vp=load([mdldir,'crust1.vp']);
vs=load([mdldir,'crust1.vs']);
rho=load([mdldir,'crust1.rho']);
bnd=load([mdldir,'crust1.bnds']);
vp1=zeros(nd,nla,nlo);
vs1=zeros(nd,nla,nlo);
rho1=zeros(nd,nla,nlo);
bnd1=zeros(nd,nla,nlo);
nline=1;
for j = 1:nla
   for i = 1:nlo
      for k = 1:nd
         vp1(k,j,i)=vp(nline,k);
         vs1(k,j,i)=vs(nline,k);
         rho1(k,j,i)=rho(nline,k);
         bnd1(k,j,i)=bnd(nline,k);
      end
      nline=nline+1;
   end
end
ilat = int16(floor(90.-lat))+1;
ilon = int16(floor(180.+lon))+1;
vptmp=vp1(:,ilat,ilon);
vstmp=vs1(:,ilat,ilon);
rhotmp=rho1(:,ilat,ilon);
ztmp=bnd1(:,ilat,ilon);
% remove the first two layers since they are water and ice, start with
% the third layer
% sediments is from layer 3 to layer 5, and crust is from layer 6 to 8
% remember converts the crustal depth to thickness (total thickness), such
% that the depth starts at 0 km
if if_topo == 0
    model = [-(ztmp(3:end)-ztmp(1)) vptmp(3:end) vstmp(3:end) rhotmp(3:end)];
else
    model = [-ztmp(3:end) vptmp(3:end) vstmp(3:end) rhotmp(3:end)];
end
% remove the layer with 0 velcoity or density (non-existing layer) as well
% as 0 thickness
thick = diff(model(:,1));
% deal with the upper mantle, assuming its a half space
thick(end+1) = 999;
knock_out = model(:,2) == 0 |  model(:,3) == 0 | model(:,4) == 0 | thick == 0;
% find the number sedimentary layers that exist
nsedi = sum(vptmp(3:5)>0 & thick(1:3) >0);
model(knock_out,:) = [];
% set the Moho depth to H
if H ~= 0
    model(end,1) = H;
end
varargout{1} = model;
varargout{2} = nsedi;
end
