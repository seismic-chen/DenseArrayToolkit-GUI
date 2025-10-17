function MigResult = leastSquaresMig3D(gather, gridStruct, param)

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



% If you need to compute great-circle distances, do so here,
% or if it's a simple local coordinate system, you can handle that as well.

% Unpack some velocity model fields
x    = gridStruct.x;
y    = gridStruct.y;
z    = gridStruct.z;
dx   = gridStruct.dx;
dy   = gridStruct.dy;
dz   = gridStruct.dz;
vp   = gridStruct.vp;
vs   = gridStruct.vs;
nx   = gridStruct.nx;
ny   = gridStruct.ny;
nz   = gridStruct.nz;
xpad = abs(x(1));  % Horizontal padding if relevant
zmax = max(z);

% If you want to restrict the imaging domain according to zmax:
zmaxSamples = floor(param.zmax / dz);
if zmaxSamples < nz
    % You can truncate the velocity arrays or just keep them as is
    % For demonstration, we won't forcibly truncate them in this example
end

%% ------------------------------------------------------------------------
%  (2) Extract trace data from gather (RF data and geometry)
% -------------------------------------------------------------------------

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

% Combine all itr traces into a matrix [Nt x Ntrace]
itrMat = cell2mat(cellfun(@(rf) rf.itr, itrAll, 'UniformOutput', false)); 

% Time axis from the first valid record
timeAxis = gval(1).RF.ittime;  % or gather(find(validMask,1)).RF.ittime
dt_samp  = gval(1).TimeAxis.dt_resample;
nt       = length(timeAxis);

%% ------------------------------------------------------------------------
%  (5) Forward wavefield modeling (source, shift, etc.)
% -------------------------------------------------------------------------
% Gather average ray param & baz if needed
raypAll = [gval.TravelInfo];
avgRayp  = mean([raypAll.rayParam]); 
% avgBaz   = mean([raypAll.baz]);
vpBottom   = mean(vp(end,:));  % example: near bottom row or top row, depending

take_off = asind(avgRayp*vpBottom/(6371-zmax)); % incident angle
back_azimuth = mean([raypAll.baz]);

timelag = calculate_timeshift(take_off,back_azimuth,x,y,vpBottom); % for plane wave time delay
%% ------------------------------------------------------------------------
%  source time function
% -------------------------------------------------------------------------
srctmp = mean(itrMat,2);

[win] = filt_win(srctmp,timeAxis,-2,1,-1);
src=srctmp.*win;
src_func = src./max(src); 

% taper
[win] = filt_win(itrMat,timeAxis,-2,20,0.5);
win = win*ones(1,size(itrMat,2));
itr=itrMat.*win;


%% ------------------------------------------------------------------------
%  (6) shifting reveiver function
% -------------------------------------------------------------------------
vel = repmat(vp,1,1,ny);
vel_s = repmat(vs,1,1,ny);
bc = 1;

[dshift] = rflsm_shift_rfs_3D(itr,timeAxis,vel,vel_s,timelag,src_func,nt,dt_samp,bc,param);


%% ------------------------------------------------------------------------
%  (6) Migration
% -------------------------------------------------------------------------
[mig, ~] = rflsm_migration_3D(dshift,timelag, vel, vel_s, src_func, nt, dt_samp,bc, param);

%% ------------------------------------------------------------------------
%  (7) LSM
% -------------------------------------------------------------------------
% We can refine param if needed
param.itermax = param.itermax;  % just reaffirming usage
param.mu      = param.mu;

[migls, ~] = rflsm_lsm_3D(dshift,timelag, vel, vel_s, src_func, nt, dt_samp,bc, param);

%% ------------------------------------------------------------------------
%  (8) Plot migration results
if param.plotMig
    % -------------------------------------------------------------------------
    figure('Name','Migration Results','Color','w','Position',[100, 100, 800, 800]);
    % Normal Migration
    subplot(2,1,1);
    imagesc(x, z, mig(:,:,5) ./ max(abs(mig(:)))); hold on;
    axis([x(1) x(end) 0 100]);
    xlabel('Distance (km)'); ylabel('Depth (km)');
    title('Migration Image');
    set(gca,'FontSize',14);
    caxis([-0.2 0.2]); colormap(seismic(3)); colorbar;

    % LSM
    subplot(2,1,2);
    imagesc(x, z, migls(:,:,10) ./ max(abs(migls(:)))); hold on;
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
MigResult.y     = y;
MigResult.z     = z;
MigResult.mig   = mig;
MigResult.migls = migls;
MigResult.param = param;

end

