function [DataStruct,noresult] = deconv_single(DataStruct, param)
% DECONV  Perform receiver function (RF) deconvolution on [T, R, Z] waveforms,
%         storing the results in the DataStruct.RF field.
%
% Usage:
%   DataStruct = deconv(DataStruct, param)
%
% Inputs:
%   DataStruct : struct array with fields:
%       Waveforms.dataProcessed -> [Nt x 3] = [T, R, Z] (from preprocessing)
%       TimeAxis.t_resample     -> time axis matching dataProcessed
%       ProcHistory             -> cell array for logging (recommended)
%
%   param : struct of deconvolution parameters.  Default fields include:
%       .gauss       = 2.5     (Gaussian width)
%       .waterlevel  = 0.01    (water-level)
%       .itmax       = 100     (max iterations for iterative decon)
%       .minderr     = 1e-5    (minimum residual threshold)
%       .phaseshift  = 5       (shift in seconds before P-wave for alignment)
%       .verbose     = false   (verbosity within decon routines)
%       .radonfilter = false   (whether to use radon-filtered data)
%
% Output:
%   DataStruct : the same input array but each valid record now contains
%       DataStruct.RF with fields:
%         .wlr      = Water-level RF trace
%         .wlrms    = RMS of water-level decon residual
%         .nwl      = Additional outputs from makeRFwater_ammon (optional)
%         .wltime   = Time axis for water-level RF
%         .itr      = Iterative decon RF trace
%         .itrms    = Final RMS or iteration error
%         .ittime   = Time axis for iterative decon
%
%   Logs are appended to DataStruct.ProcHistory for each step.
%
% Example:
%   config.DeconvParam.gauss = 2.5;
%   config.DeconvParam.itmax = 200;
%   DataStruct = deconv(DataStruct, config.DeconvParam);
%
% Author: Yunfeng Chen (Refined by ChatGPT)
% Date  :  Feb. 3, 2025

%% 1) Set default parameters if missing
if nargin < 2, param = struct(); end
if ~isfield(param, 'gauss'),       param.gauss       = 2.5;   end
if ~isfield(param, 'waterlevel'),  param.waterlevel  = 0.01;  end
if ~isfield(param, 'itmax'),       param.itmax       = 100;   end
if ~isfield(param, 'minderr'),     param.minderr     = 1e-5;  end
if ~isfield(param, 'phaseshift'),  param.phaseshift  = 5;     end
if ~isfield(param, 'verbose'),     param.verbose     = false; end
if ~isfield(param, 'radonfilter'), param.radonfilter = false; end

% noresult will store flag of records that fail decon or produce no result
noresult = false;

%% 2.1 Check for processed waveforms
waveKey = 'dataProcessed';  % default: use the preprocessed waveforms
if param.radonfilter
    % If radonfilter = true, we prefer "dataRadonFiltered"
    waveKey = 'dataRadonFiltered';
end

if ~isfield(DataStruct.Waveforms, waveKey) || ...
        isempty(DataStruct.Waveforms.(waveKey))
    warnMsg = sprintf('[Deconv] No %s waveforms for trace -> skip.', waveKey);
    warning(warnMsg);
    DataStruct = appendHistory(DataStruct, warnMsg);
    noresult=true;
    return;
end

seisPRZ = DataStruct.Waveforms.(waveKey);  % e.g., [Nt x 3] = [T, R, Z]
if size(seisPRZ,2) < 3
    warnMsg = sprintf('[Deconv] %s has <3 comps for trace -> skip.', waveKey);
    warning(warnMsg);
    DataStruct = appendHistory(DataStruct, warnMsg);
    noresult=true;
    return;
end

%% 2.2 Check time axis
if ~isfield(DataStruct.TimeAxis, 't_resample') || ...
        isempty(DataStruct.TimeAxis.t_resample)
    warnMsg = sprintf('[Deconv] t_resample missing for trace -> skip.');
    warning(warnMsg);
    DataStruct = appendHistory(DataStruct, warnMsg);
    noresult=true;
    return;
end

t = DataStruct.TimeAxis.t_resample;
if length(t) ~= size(seisPRZ,1)
    warnMsg = sprintf('[Deconv] length(t)=%d but size(seisPRZ,1)=%d mismatch -> skip trace', ...
        length(t), size(seisPRZ,1));
    warning(warnMsg);
    DataStruct = appendHistory(DataStruct, warnMsg);
    noresult=true;
    return;
end

%% 2.3 Extract R/Z components (assuming T=1, R=2, Z=3)
R  = seisPRZ(:,2);
Z  = seisPRZ(:,3);
dt = t(2) - t(1);
nt = length(t);

% Make sure DataStruct.RF is at least an empty struct to fill results
if ~isfield(DataStruct, 'RF'), DataStruct.RF = struct(); end

%% 2.4 Water-level deconvolution
try
    % (The function signature of makeRFwater_ammon may vary in your codebase)
    [wlr, wlrms, nwl] = makeRFwater_ammon( ...
        R, Z, param.phaseshift, dt, nt, param.waterlevel, param.gauss, param.verbose);

    % Store in DataStruct
    DataStruct.RF.wlr    = wlr(:);
    DataStruct.RF.wlrms  = wlrms;
    DataStruct.RF.nwl    = nwl;
    DataStruct.RF.wltime = (dt*(0:nt-1) - param.phaseshift)';

    % Log message
    logMsg = sprintf('[Deconv] Water-level decon OK (gauss=%.2f, wlevel=%.3f)', ...
        param.gauss, param.waterlevel);
    DataStruct = appendHistory(DataStruct, logMsg);

catch ME
    warnMsg = sprintf('[Deconv] Trace: Water-level decon failed: %s', ME.message);
    warning(warnMsg);
    DataStruct = appendHistory(DataStruct, warnMsg);
    noresult=true;
    return;
end

%% 2.5 Iterative deconvolution
try
    % (Again, adapt to your actual function signature for iterative decon)
    [itr, itrms] = makeRFitdecon_la_norm( ...
        R, Z, dt, nt, param.phaseshift, param.gauss, param.itmax, param.minderr);

    if isempty(itrms)
        % If itrms is empty, no valid solution
        warnMsg = sprintf('[Deconv] Trace: empty itrms -> no iterative RF.');
        warning(warnMsg);
        DataStruct = appendHistory(DataStruct, warnMsg);
        noresult=true;
        return;
    else
        % Store iterative result
        DataStruct.RF.itr    = itr(:);
        DataStruct.RF.itrms  = itrms(end);
        DataStruct.RF.ittime = (dt*(0:nt-1) - param.phaseshift)';

        logMsg = sprintf('[Deconv] Iterative decon OK (gauss=%.2f, itmax=%d, minderr=%.1e)', ...
            param.gauss, param.itmax, param.minderr);
        DataStruct = appendHistory(DataStruct, logMsg);
    end
catch ME
    warnMsg = sprintf('[Deconv] Trace: Iterative decon failed: %s', ME.message);
    warning(warnMsg);
    DataStruct = appendHistory(DataStruct, warnMsg);
    noresult=true;
    return;
end
end

%% Helper function: appendHistory
function DataStruct = appendHistory(DataStruct, msg)
% APPENDHISTORY  Append a message to DataStruct.ProcHistory.
%                If it does not exist, create it as a cell array.

if ~isfield(DataStruct, 'ProcHistory') || isempty(DataStruct.ProcHistory)
    DataStruct.ProcHistory = {msg};
else
    DataStruct.ProcHistory{end+1} = msg;
end
end