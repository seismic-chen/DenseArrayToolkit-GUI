function ccpResult = CCPCommonEventGather(gather, velocityModel, gridStruct)
% CCPCommonEventGather  Perform Common Conversion Point stacking (CCP) on seismic data.
%
% Usage:
%   ccpResult = CCPCommonEventGather(gather, velocityModel, gridStruct)
%
% Inputs:
%   gather         : Struct array of seismic traces for one event or gather.
%                    Each element includes fields:
%                       .RF.itr            - Iterative decon RF [Nt x 1]
%                       .TravelInfo.rayParam, .TravelInfo.baz - Ray parameters
%                       .TimeAxis.t_resample, .TimeAxis.dt_resample - time axis
%   velocityModel  : Struct describing velocity model with fields:
%                       .x, .z, .nx, .nz, .dx, .dz, .vp, .vs, ...
%   gridStruct  : Struct with profile-related info:

%
% Outputs:
%   ccpResult : Struct containing CCP results with fields:
%       .x, .z    - Horizontal and depth axes
%       .img      - CCP image
%       .param    - A copy of the param struct for reference
%
% Author: Yunfeng Chen
% Date: Feb. 18, 2025

%% Input Validation and Parameter Defaults
% Validate input data structure and parameters
% --------------------------------------------------
% Validate seismic gather structure
% - Check if gather is non-empty struct array
% - Each gather element should contain RF data and travel time information
if isempty(gather) || ~isstruct(gather)
    error('CCPCommonEventGather:InvalidGather', 'Gather must be a non-empty struct array.');
end

% Validate velocityModel fields (currently disabled)
% requiredVMfields = {'x', 'z', 'vp', 'vs', 'dx', 'dz', 'nx', 'nz'};
% Uncomment and complete this validation if needed
% for f = requiredVMfields
%     if ~isfield(velocityModel, f{1})
%         error('CCPCommonEventGather:MissingVelocityModelField', 'Field "%s" is missing in velocityModel.', f{1});
%     end
% end

% Display event info if available
if isfield(gather(1), 'EventInfo') && isfield(gather(1).EventInfo, 'evid')
    disp(['Processing event: ', gather(1).EventInfo.evid]);
end

%% Unpack Profile and Velocity Model Information
% Extract velocity model parameters and grid configuration
% --------------------------------------------------
% velocityModel contains either 1D or 3D velocity information
% gridStruct defines the imaging grid parameters
if strcmp(velocityModel.ModelType ,'1D')
    vp   = velocityModel.vp(:, 1);  % Assuming vp is 1D
    vs   = velocityModel.vs(:, 1);  % Assuming vs is 1D
    z = velocityModel.z;
else
    [z, r, vp, vs, ~, ~] = ak135('cont');
end
dz = gridStruct.dz;
zmax = max(gridStruct.z);
zout = 0:dz:zmax;
%% Calculate Conversion Points and Ray Tracing
% --------------------------------------------------
% Key Steps:
% 1. Extract receiver function (RF) data and time axes
% 2. Calculate ray parameters from travel time information
% 3. Perform ray tracing to find conversion points
% 4. Apply 3D velocity corrections if needed
nrf = length(gather);

% Extract RFs and times from gather
rfsAll = cellfun(@(rf) rf.itr, {gather.RF}, 'UniformOutput', false);
timeAll = cellfun(@(rf) rf.ittime, {gather.RF}, 'UniformOutput', false);

% Extract ray parameters and back-azimuths from gather
raypAll = cellfun(@(ti) ti.rayParam / 6371, {gather.TravelInfo}, 'UniformOutput', false);
raypAll = cell2mat(raypAll);
bazAll = cellfun(@(ti) ti.baz, {gather.TravelInfo}, 'UniformOutput', false);
bazAll = cell2mat(bazAll);

% Extract station coordinates from gather
latAll = cellfun(@(si) si.stla, {gather.StationInfo}, 'UniformOutput', false);
latAll = cell2mat(latAll);
lonAll = cellfun(@(si) si.stlo, {gather.StationInfo}, 'UniformOutput', false);
lonAll = cell2mat(lonAll);

% Perform ray tracing for CCP points
model_type = 'flat';
disp('Ray tracing started');
tic;
[cp, ~ ,MidPoints] = rf_ccp(raypAll, bazAll, gridStruct.dz, zmax, z, vp, vs, latAll, lonAll, model_type);
toc;
disp('Ray tracing completed');

% make time correction if a regional 3D velocity model is available
if strcmp(velocityModel.ModelType ,'3D')
    Fvp = velocityModel.vp;
    Fvs = velocityModel.vs;
    % correct for heterogeneity
    RayDepths = (1*dz:dz:zmax)';
    [TimeCorrections, ~, ~] = correct_RFs(MidPoints, RayDepths, Fvp, Fvs, z, vp, vs);
else
    TimeCorrections = zeros(length(zout),length(rfsAll)); 
end

% Time-to-Depth Conversion
% --------------------------------------------------
% Convert RFs from time domain to depth domain using:
% - Ray parameters (raypAll)
% - Velocity model (vp, vs)
% - Time corrections from 3D velocity model (if applicable)
% This converts the time-based RF measurements to depth coordinates
% matching the imaging grid
disp('Time-to-depth conversion started');
tic;
[~, rfsAll_depth, ~] = rf_migrate(timeAll, rfsAll, raypAll, gridStruct.dz, zmax, z, vp, vs, TimeCorrections);
toc;
disp('Time-to-depth conversion completed');


% Attach RF amplitudes to CCP points
for k = 1:nrf
    cp(k).amp = rfsAll_depth{k};
end

% Project CCP points to profile
for k = 1:nrf
    slon_tmp = [cp(k).lonb];
    slat_tmp = [cp(k).latb];
    [rx, ry] = latlonToProjectedCoords(slon_tmp, slat_tmp, gridStruct);
    cp(k).rx = rx;
    cp(k).ry = ry;
end

%% CCP Stacking Process
% --------------------------------------------------
% Core algorithm steps:
% 1. Initialize imaging grid based on gridStruct parameters
% 2. Bin RF amplitudes into spatial grid cells
% 3. Stack (sum) amplitudes in each cell
% 4. Normalize by sample count in each cell
% 
% Handles both 1D and 3D velocity models differently:
% - 1D: Simple 2D (distance-depth) stacking
% - 3D: Full 3D spatial binning with progress tracking
switch velocityModel.ModelType
    % project to profile for 2D imaging
    case '1D'
        % Set grid sizes (distance direction: dx, depth direction: dz)
        [X, Z] = meshgrid(gridStruct.x, gridStruct.z);
        nx = length(gridStruct.x);
        nz = length(gridStruct.z);

        % Initialize stacking result matrices
        V = zeros(nz, nx); % Store accumulated amplitude values
        count = zeros(nz, nx); % Store sample count for each grid
        % Perform parallel processing for speed
        parfor i = 1:nz
            for j = 1:nx
                for n = 1:length(cp)
                    xi = cp(n).rx;
                    zi = cp(n).zpos;
                    vi = cp(n).amp;
                    keep = xi >= gridStruct.x(j) - 2 * gridStruct.dx & xi <= gridStruct.x(j) + 2 * gridStruct.dx & ...
                        zi >= gridStruct.z(i) - 2 * gridStruct.dz & zi <= gridStruct.z(i) + 2 * gridStruct.dz;
                    V(i, j) = V(i, j) + sum(vi(keep));
                    count(i, j) = count(i, j) + sum(keep);
                end
            end
        end

        % Generate CCP Image
        V = V ./ max(count, 1);  % 避免除以0

        %% Plot CCP Image
        try
            load roma;
            cmap = flipud(roma);
        catch
            cmap = parula;
        end

        % figure;
        % set(gcf, 'Position', [100 100 800 400], 'color', 'w');
        % imagesc(gridStruct.x, gridStruct.z, V);
        % caxis([-0.1 0.1]);
        % colormap(cmap);
        % colorbar;
        % xlabel('Distance (km)');
        % ylabel('Depth (km)');
        % title('CCP image');
        % set(gca, 'fontsize', 14);

        %% Save CCP results
        ccpResult = struct('X', X, 'Z', Z, 'img', V, 'count', count);
    case '3D'
        nx = gridStruct.nx;
        ny = gridStruct.ny;
        nz = gridStruct.nz;
        % Set grid sizes (distance direction: dx, depth direction: dz)
        [X, Y, Z] = meshgrid(gridStruct.x, gridStruct.y, gridStruct.z);

        % Initialize stacking result matrices
        V = zeros(ny, nx, nz); % Store accumulated amplitude values
        count = zeros(ny, nx, nz); % Store sample count for each grid
        xmin = min(gridStruct.x);
        ymin = min(gridStruct.y);
        dx = gridStruct.dx;
        dy = gridStruct.dy;

        for n=1:length(cp)
            if mod(n,10) == 0
                disp(['Binning ',num2str(n),'/',num2str(length(cp)),' traces']);
            end
            xx=cp(n).rx;
            yy=cp(n).ry;
            zz=cp(n).zpos;
            
            % 3D stacking
            for k=1:length(zz)
                i=floor((yy(k)-ymin)/dy)+1;
                j=floor((xx(k)-xmin)/dx)+1;
        % Calculate z-layer index
        m = floor((zz(k) - gridStruct.z(1)) / gridStruct.dz) + 1;
        % Validate array indices
                if i>=1 && i<=ny && j>=1 && j<=nx && m>=1 && m<=nz
                    amp = cp(n).amp(k);
                    if ~isnan(amp)
                        count(i,j,m) = count(i,j,m) + 1;
                        V(i,j,m) = V(i,j,m) + amp;
                    end
                end
            end

        end
%         V = V./max(count,1);  % Disabled alternative normalization

        ccpResult = struct('X', X, 'Y', Y, 'Z', Z, 'img', V, 'count', count);
end
end
