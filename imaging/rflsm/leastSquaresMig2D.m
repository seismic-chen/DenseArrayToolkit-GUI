function MigResult = leastSquaresMig2D(gather, gridStruct, param)
% LEASTSQUARESMIG  Perform least squares migration (LSM) on seismic data.
%
% Usage:
%   MigResult = leastSquaresMig2D(gather, gridStruct, param)
%
% Inputs:
%   gather         : A struct array of seismic traces for one event or gather.
%                    Each element typically includes:
%                       .RF.itr      - Iterative decon RF [Nt x 1]
%                       .TravelInfo.rayParam, .baz - Ray parameter, back-azimuth
%                       .TimeAxis.t_resample, .TimeAxis.dt_resample - time axis
%
%
%   gridStruct  : Struct with grid-related info, e.g.:
%
%   param          : Migration and reconstruction parameters, for example:
%       .ssa        (bool)   - if true, use SSA-based reconstruction
%       .dz         (double) - depth sampling for migration (km)
%       .zmax       (double) - maximum depth for imaging (km)
%       .binning.dx (double) - bin size in horizontal distance
%       .plotBinned (bool)   - if true, show binned data wiggle plot
%       .itermax    (int)    - max iteration for LSM
%       .mu         (double) - damping or regularization factor
%       .xpad       (double) - horizontal padding (for advanced modeling)
%       .plotMig    (bool)   - if true, plot migration results
%
% Outputs:
%   MigResult : Struct containing migration results. Common fields:
%       .x, .z    - horizontal and depth axes
%       .mig      - migrated image (normal migration)
%       .migls    - LSM result
%       .param    - a copy of the param struct for reference
%
% Example:
%   config.MigParam.zmax = 100;
%   MigResult = leastSquaresMig(gather, profile, config.MigParam);
%
% Author:  Yunfeng Chen (Refined by ChatGPT)
% Date   :  Jan. 25, 2025

%% ------------------------------------------------------------------------
%  (0) Check input structures & fill default parameters
% -------------------------------------------------------------------------
if nargin < 3,  param = struct(); end

% Ensure gather has at least one valid entry
if isempty(gather) || ~isstruct(gather)
    error('leastSquaresMig:InvalidGather', ...
          'Gather is empty or not a struct array.');
end

% Fill param defaults
if ~isfield(param,'ssa'),        param.ssa = false; end
if ~isfield(param,'dz'),         param.dz = 1;      end
if ~isfield(param,'zmax'),       param.zmax = 800;  end
if ~isfield(param,'plotBinned'), param.plotBinned = false; end
if ~isfield(param,'binning') || ~isfield(param.binning,'dx')
    % If binning not given, default to gridStruct.dx
    param.binning = gridStruct;
end
if ~isfield(param,'itermax'), param.itermax = 20; end
if ~isfield(param,'mu'),      param.mu = 0.1;     end
if ~isfield(param,'xpad'),    param.xpad = 0;     end
if ~isfield(param,'plotMig'), param.plotMig = false;    end
% Display some gather info if present
if isfield(gather(1), 'EventInfo') && isfield(gather(1).EventInfo, 'evid')
    disp(['[leastSquaresMig] Processing event: ' gather(1).EventInfo.evid]);
end

%% ------------------------------------------------------------------------
%  (1) Unpack profile & velocity model info
% -------------------------------------------------------------------------
% Example usage of profileStruct: (lon1, lat1) to (lon2, lat2)
lon1 = gridStruct.principalAxisLatLon(1, 1);
lat1 = gridStruct.principalAxisLatLon(1, 2);
lon2 = gridStruct.principalAxisLatLon(end, 1);
lat2 = gridStruct.principalAxisLatLon(end, 2);

evla = gather(1).EventInfo.evla;
evlo = gather(1).EventInfo.evlo;

% check distance to end points of the profile
[dist1,az1] = distance(evla,evlo,lat1,lon1);
[dist2,az2] = distance(evla,evlo,lat2,lon2);

% incidence direction is 1 if wavefield is propagating towards the
% positive direction of the profile and -1 otherwise
if dist2 > dist1
    incidence_direction = 1;
else
    incidence_direction = -1;
end

% project to profile
stla = cellfun(@(stationinfo) stationinfo.stla, {gather.StationInfo}, 'UniformOutput', false);
stlo = cellfun(@(stationinfo) stationinfo.stlo, {gather.StationInfo}, 'UniformOutput', false);
stlo = cell2mat(stlo)';
stla = cell2mat(stla)';
[rx, ry] = latlonToProjectedCoords(stlo, stla, gridStruct);
for n=1:length(gather)
    gather(n).RF.rx = rx(n);
end

% If you need to compute great-circle distances, do so here,
% or if it's a simple local coordinate system, you can handle that as well.

% Unpack some velocity model fields
x    = gridStruct.x;
z    = gridStruct.z;
dx   = gridStruct.dx;
dz   = gridStruct.dz;
vp   = gridStruct.vp;
vs   = gridStruct.vs;
nx   = gridStruct.nx;
nz   = gridStruct.nz;
xpad = abs(x(1));  % Horizontal padding if relevant

% If you want to restrict the imaging domain according to zmax:
zmaxSamples = floor(param.zmax / dz);
if zmaxSamples < nz
    % You can truncate the velocity arrays or just keep them as is
    % For demonstration, we won't forcibly truncate them in this example
end

%% ------------------------------------------------------------------------
%  (2) Extract trace data from gather (RF data and geometry)
% -------------------------------------------------------------------------
%  gather(i).RF.itr => iterative decon RF [Nt x 1]
%  gather(i).RF.rx  => receiver x-loc if defined (or station offset in x)
%  gather(i).TravelInfo.rayParam, .baz => used if needed
%  gather(i).TimeAxis.t_resample, dt_resample => time axis

% Check which gather entries have .RF.itr
itrAll = {gather.RF};  % cell array of RF structures
validMask = cellfun(@(r) ~isempty(r) && isfield(r,'itr') && ~isempty(r.itr), itrAll);

if ~any(validMask)
    warning('[leastSquaresMig] No valid .RF.itr found in gather. Returning empty result.');
    MigResult = struct();
    return;
end

% Keep only valid elements
itrAll  = itrAll(validMask);
gval    = gather(validMask);

% Extract Rx positions (assuming gather(i).RF.rx is the horizontal location)
rxCell  = cellfun(@(rf) rf.rx, itrAll, 'UniformOutput', false);
rx      = cell2mat(rxCell);           % [1 x Ntrace]

% Combine all itr traces into a matrix [Nt x Ntrace]
itrMat = cell2mat(cellfun(@(rf) rf.itr, itrAll, 'UniformOutput', false));

% Time axis from the first valid record
timeAxis = gval(1).RF.ittime;  % or gather(find(validMask,1)).RF.ittime
dt_samp  = gval(1).TimeAxis.dt_resample;
nt       = length(timeAxis);

%% ------------------------------------------------------------------------
%  (3) Binning the traces along the x-axis
% -------------------------------------------------------------------------
% Binning step: group traces into horizontal bins of width param.binning.dx
dBinned = doBinning(itrMat, rx, x, param.binning.dx);
% dBinned => [Nt x nx], each column is the average of all traces that fall
%            into that bin's x-range

% Optional: Plot binned data if requested
if param.plotBinned
    figure('Name','Binned Data','Color','w');
    wigb(dBinned, 2, x, timeAxis);
    ylim([-5 30]);  % Adjust as necessary
    xlabel('Distance (km)'); ylabel('Time (s)');
    title('Binned RF Data');
end

%% ------------------------------------------------------------------------
%  (4) Optional: SSA or advanced reconstruction
% -------------------------------------------------------------------------
if param.ssa
    % For example, skip a padding region if xpad>0
    % i1, i2 define subrange in x?
    i1 = floor(xpad/dx) + 1;
    i2 = nx - floor(xpad/dx);
    i2 = min(i2, size(dBinned,2));  % just to be safe

    if i1 < i2
        din  = dBinned(:, i1:i2);
        % Call your custom function rflsm_ssa() or other advanced method:
        dSSA = rflsm_ssa(din, dt_samp, param);
        % Place reconstructed data back
        dBinned(:, i1:i2) = dSSA;
        disp('[leastSquaresMig] SSA-based reconstruction done.');
    end
end

%% ------------------------------------------------------------------------
%  (5) Forward wavefield modeling (source, shift, etc.)
% -------------------------------------------------------------------------
% Gather average ray param & baz if needed
raypAll = [gval.TravelInfo];
avgRayp  = mean([raypAll.rayParam]) / 6371;  % e.g., if rayParam is in s/deg, or s/radius
% avgBaz   = mean([raypAll.baz]);
vpBottom   = mean(vp(end,:));  % example: near bottom row or top row, depending

% Create or define a source wavelet for LSM
[src, pos, tshift] = rflsm_create_src(dt_samp, nt, avgRayp, incidence_direction, vpBottom, param);
% Scale source if needed
src = src * max(mean(itrMat, 2));

% Shift data if needed by travel-time differences
dshift = rflsm_shift_rfs(dBinned, timeAxis, vp, vs, src, pos, tshift, param);

%% ------------------------------------------------------------------------
%  (6) Migration
% -------------------------------------------------------------------------
[mig, ~] = rflsm_migration(dshift, timeAxis, vp, vs, src, pos, tshift, param);

%% ------------------------------------------------------------------------
%  (7) LSM
% -------------------------------------------------------------------------
% We can refine param if needed
param.itermax = param.itermax;  % just reaffirming usage
param.mu      = param.mu;

[migls, ~] = rflsm_lsm(dshift, timeAxis, vp, vs, src, pos, tshift, param);

%% ------------------------------------------------------------------------
%  (8) Plot migration results
if param.plotMig
    % -------------------------------------------------------------------------
    figure('Name','Migration Results','Color','w','Position',[100, 100, 800, 800]);

    % Normal Migration
    subplot(2,1,1);
    imagesc(x, z, mig ./ max(abs(mig(:)))); hold on;
    axis([x(1) x(end) 0 100]);
    xlabel('Distance (km)'); ylabel('Depth (km)');
    title('Migration Image');
    set(gca,'FontSize',14);
    caxis([-0.2 0.2]); colormap(seismic(3)); colorbar;

    % LSM
    subplot(2,1,2);
    imagesc(x, z, migls ./ max(abs(migls(:)))); hold on;
    axis([x(1) x(end) 0 100]);
    xlabel('Distance (km)'); ylabel('Depth (km)');
    title('Least Squares Migration (LSM)');
    set(gca,'FontSize',14);
    caxis([-0.2 0.2]); colorbar;
end

%% ------------------------------------------------------------------------
%  (9) Prepare output structure
% -------------------------------------------------------------------------
MigResult.x     = x;
MigResult.z     = z;
MigResult.mig   = mig;
MigResult.migls = migls;
MigResult.param = param;

end

%% ========================================================================
%  SUBFUNCTIONS
%% ========================================================================

function dBinned = doBinning(itrMat, rx, xgrid, dx)
% DOBINNING  Bin trace data (itrMat) along x, using bins centered at xgrid.
%
% Inputs:
%   itrMat  : [Nt x Ntrace], each column is a trace
%   rx      : [1 x Ntrace], horizontal location of each trace
%   xgrid   : array of bin centers [1 x nx]
%   dx      : bin width (km)
%
% Output:
%   dBinned : [Nt x nx], binned (averaged) data

[nt, nTr] = size(itrMat);
nx = length(xgrid);
dBinned = zeros(nt, nx);

for i = 1:nx
    xLeft  = xgrid(i) - dx/2;
    xRight = xgrid(i) + dx/2;
    inBin  = (rx >= xLeft) & (rx < xRight);
    if any(inBin)
        % Average all traces that fall into this bin
        dBinned(:, i) = mean(itrMat(:, inBin), 2);
    end
end

end

