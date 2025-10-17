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

param.paramMig.Ti = timeAxis;
param.paramMig.dt = dt_samp;
param.paramMig.nt = nt;
param.paramMig.gauss = param.gauss;
param.paramMig.phaseshift = param.phaseshift;
%% ------------------------------------------------------------------------
%  (5) Forward wavefield modeling (source, shift, etc.)
% -------------------------------------------------------------------------
% Gather average ray param & baz if needed
raypAll = [gval.TravelInfo];
avgRayp  = mean([raypAll.rayParam]); 
% avgBaz   = mean([raypAll.baz]);
vp = param.paramMig.vp;
vpBottom   = mean(vp(end,:));  % example: near bottom row or top row, depending

take_off = asind(avgRayp*vpBottom/(6371-max(gridStruct.z))); % incident angle
back_azimuth = mean([raypAll.baz]);

disp(['Event: Baz ', num2str(back_azimuth),', Inc: ',num2str(take_off)])
disp('---------------------------------------------------------------')
%% shift RF
stla = cell2mat(cellfun(@(sta) sta.stla, {gather.StationInfo}, 'UniformOutput', false));    
stlo = cell2mat(cellfun(@(sta) sta.stlo, {gather.StationInfo}, 'UniformOutput', false));    

[rx, ry] = latlonToProjectedCoords(stlo, stla, gridStruct);
x = gridStruct.x;
y = gridStruct.y;
xo = gridStruct.XInOriginalCoord(:);
yo = gridStruct.YInOriginalCoord(:);
[rfshift,src_func,mask] = shiftRFs(itrMat,take_off,back_azimuth,xo,yo,x,y,rx,ry,param);

%% migration
save_wavefield = 0;
[mig,pre_rfm] = runMigration(rfshift,take_off,back_azimuth,src_func,save_wavefield,gridStruct,param);

%% LSM
itermax = 15;
save_wavefield = 0;
[lsmig,pre_rflsm] = runLSM(rfshift,take_off,back_azimuth,src_func,save_wavefield,itermax,gridStruct,param);

%% ------------------------------------------------------------------------
%  (8) Prepare output structure
% -------------------------------------------------------------------------
MigResult.mig   = mig;
MigResult.migls = lsmig;
x = param.paramMig.x;
y = param.paramMig.y;
z = param.paramMig.z;
[X,Y,Z] = meshgrid(x,y,z);
MigResult.X = X;
MigResult.Y = Y;
MigResult.Z = Z;
%% ------------------------------------------------------------------------
%  (9) Plot migration results

% plotResults(mig,pre_rfm,lsmig,pre_rflsm,param,mask);


end

