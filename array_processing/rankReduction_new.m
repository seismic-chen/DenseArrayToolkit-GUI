function [gather, d1_otg] = rankReduction_new(gather, gridStruct, param)
% RANKREDUCTION - Do rank reduction (DRR-OTG) on gather data in 3D (time x x-dist x y-dist)
%
% Usage:
%   [gather, d1_otg] = rankReduction(gather, param)
%
% Inputs:
%   gather - struct array of seismic records, each element typically has:
%       .StationInfo.stlo, .StationInfo.stla
%       .RF.itr         => [Nt x 1], the RF trace
%       .RF.ittime      => [Nt x 1], time axis
%       .TimeAxis.dt_resample => sampling interval (or .dt?)
%   param  - parameter struct with fields:
%       .lonmin, .latmin : reference lat/lon for coordinate transform
%       .nx, .ny         : grid dimensions
%       .rank            : rank for DRR
%       .K               : damping factor (from your DRR method)
%       .niter           : iteration count
%       .eps             : small regularization or tolerance
%       .verb            : verbose mode (bool)
%       .mode            : DRR mode or algorithm control
%       .flow, .fhigh    : frequency bounds
%
% Outputs:
%   gather - updated gather, where gather(n).RF.itr is replaced by the reconstructed trace
%   d1_otg - the 3D data volume in time-x-y after DRR-OTG reconstruction (if same size as d0)
%
% Dependencies:
%   getStations, latlon2xy, drr_bin3d, drr3drecon_otg 
%
% Author: Yunfeng Chen
% Date:   Jan. 27, 2025

%% 0. Basic checks
if ~isfield(param,'nx') || ~isfield(param,'ny')
    error('rankReduction:MissingParam','param.nx and param.ny must be specified.');
end
if ~isfield(param,'flow'),   param.flow   = 0.1; end
if ~isfield(param,'fhigh'),  param.fhigh  = 1.2; end
if ~isfield(param,'rank'),   param.rank   = 10; end
if ~isfield(param,'K'),      param.K      = 5;  end
if ~isfield(param,'niter'),  param.niter  = 20; end
if ~isfield(param,'eps'),    param.eps    = 1e-3; end
if ~isfield(param,'verb'),   param.verb   = true; end
if ~isfield(param,'mode'),   param.mode   = 1; end
if ~isfield(param,'tmax'),   param.tmax   = 30; end
if ~isfield(param,'plotRankReduction'),   param.plotRankReduction   = false; end


%% 1. Get station info
stationList = getStations(gather);
% stationList should have fields .stlo, .stla
% Flatten stlo, stla to vector if needed
stlo = [stationList.stlo]';  
stla = [stationList.stla]';

% 每个地震台站在投影坐标系中的笛卡尔坐标
[rx, ry] = latlonToProjectedCoords(stlo, stla, gridStruct);

%% 3. Collect RF data into matrix d0
% gather(i).RF.itr => the trace
% Make sure gather has consistent time
if ~isfield(gather(1).RF,'ittime')
    error('rankReduction:NoTime','gather(1).RF.ittime is missing.');
end
t = gather(1).RF.ittime;
% see if we have dt
if isfield(gather(1).TimeAxis,'dt_resample')
    dt = gather(1).TimeAxis.dt_resample;
else
    warning('No dt_resample in gather(1).TimeAxis. Using default 0.01s');
    dt = 0.1;
end

% optionally cut at tmax
idxT = (t <= param.tmax);
t    = t(idxT);

% build data matrix d0  [Nt x Ntrace]
itrCell = {gather.RF};
validMask = cellfun(@(rf) isfield(rf,'itr') && ~isempty(rf.itr), itrCell);
itrCell  = itrCell(validMask);
% combine
d0 = cell2mat(cellfun(@(rf) rf.itr(idxT), itrCell,'UniformOutput', false));

%% 4. bin the data in 3D
% [d3d, x1, y1, mask] = drr_bin3d(d0, rx, ry, param.nx, param.ny, param.ox, param.oy, param.mx, param.my);
% miss_per = (length(find(mask==0))/length(t))/(param.nx*param.ny);

%% 5. DRR-OTG reconstruction
% call drr3drecon_otg
[d1_otg, d1] = drr3drecon_otg(...
    d0, rx, ry, ...
    param.nx, param.ny, ...
    param.ox, param.oy, param.mx, param.my, ...
    param.flow, param.fhigh, dt, ...
    param.rank, param.K, param.niter, param.eps, ...
    param.verb, param.mode);
%% 6. Update gather with reconstructed traces
% d1 => presumably same shape as d0 => Ntrace columns
% We assume the order of gather matches the order in which we constructed d0
% If so, we do:
validIdx = find(validMask);

if size(d1,2) ~= length(validIdx)
    warning('rankReduction:SizeMismatch', ...
       'd1 has %d columns but validIdx length = %d. Some gather traces not updated.', ...
        size(d1,2), length(validIdx));
end

for k = 1:length(validIdx)
    col = k;  % The column in d1
    gInd = validIdx(k);
    if col <= size(d1,2)
        gather(gInd).RF.itr    = d1(:, col);  % [Nt], reconstructed
        gather(gInd).RF.ittime = t;          % new time axis
    end
end

if param.plotRankReduction
    figure; imagesc(1:length(gather)*2,t,[d0 d1])
    caxis([-0.1 0.1])
    colormap(seismic(3))
end
end
