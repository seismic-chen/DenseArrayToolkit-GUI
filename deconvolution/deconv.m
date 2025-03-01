function DataStruct = deconv(DataStruct, param)
% DECONV  Perform receiver function (RF) deconvolution on [T, R, Z] waveforms,
%         storing the results in the DataStruct(n).RF field.
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
%       DataStruct(n).RF with fields:
%         .wlr      = Water-level RF trace
%         .wlrms    = RMS of water-level decon residual
%         .nwl      = Additional outputs from makeRFwater_ammon (optional)
%         .wltime   = Time axis for water-level RF
%         .itr      = Iterative decon RF trace
%         .itrms    = Final RMS or iteration error
%         .ittime   = Time axis for iterative decon
%
%   Logs are appended to DataStruct(n).ProcHistory for each step.
%
% Example:
%   config.DeconvParam.gauss = 2.5;
%   config.DeconvParam.itmax = 200;
%   DataStruct = deconv(DataStruct, config.DeconvParam);
%
% Author: Yunfeng Chen (Refined by ChatGPT)
% Date  :  Jan 9, 2025

%% 1) Set default parameters if missing
if nargin < 2, param = struct(); end
if ~isfield(param, 'gauss'),       param.gauss       = 2.5;   end
if ~isfield(param, 'waterlevel'),  param.waterlevel  = 0.01;  end
if ~isfield(param, 'itmax'),       param.itmax       = 100;   end
if ~isfield(param, 'minderr'),     param.minderr     = 1e-5;  end
if ~isfield(param, 'phaseshift'),  param.phaseshift  = 5;     end
if ~isfield(param, 'verbose'),     param.verbose     = false; end
if ~isfield(param, 'radonfilter'), param.radonfilter = false; end

% noresult will store indices of records that fail decon or produce no result
noresult = [];

disp('--- Start Deconvolution ---');
tStart = tic;  % track elapsed time

%% 2) Loop over each trace/record in DataStruct
for n = 1:length(DataStruct)
    % Show progress every 100 traces (adjust as needed)
    if mod(n, 100) == 0
        fprintf('Deconvolution on trace %d/%d\n', n, length(DataStruct));
    end

    %% 2.1 Check for processed waveforms
    waveKey = 'dataProcessed';  % default: use the preprocessed waveforms
    if param.radonfilter
        % If radonfilter = true, we prefer "dataRadonFiltered"
        waveKey = 'dataRadonFiltered';
    end

    if ~isfield(DataStruct(n).Waveforms, waveKey) || ...
            isempty(DataStruct(n).Waveforms.(waveKey))
        warnMsg = sprintf('[Deconv] No %s waveforms for trace %d -> skip.', waveKey, n);
        warning(warnMsg);
        DataStruct = appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    seisPRZ = DataStruct(n).Waveforms.(waveKey);  % e.g., [Nt x 3] = [T, R, Z]
    if size(seisPRZ,2) < 3
        warnMsg = sprintf('[Deconv] %s has <3 comps for trace %d -> skip.', waveKey, n);
        warning(warnMsg);
        DataStruct = appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    %% 2.2 Check time axis
    if ~isfield(DataStruct(n).TimeAxis, 't_resample') || ...
            isempty(DataStruct(n).TimeAxis.t_resample)
        warnMsg = sprintf('[Deconv] t_resample missing for trace %d -> skip.', n);
        warning(warnMsg);
        DataStruct = appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    t = DataStruct(n).TimeAxis.t_resample;
    if length(t) ~= size(seisPRZ,1)
        warnMsg = sprintf('[Deconv] length(t)=%d but size(seisPRZ,1)=%d mismatch -> skip trace %d.', ...
            length(t), size(seisPRZ,1), n);
        warning(warnMsg);
        DataStruct = appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    %% 2.3 Extract R/Z components (assuming T=1, R=2, Z=3)
    R  = seisPRZ(:,2);
    Z  = seisPRZ(:,3);
    dt = t(2) - t(1);
    nt = length(t);

    % Make sure DataStruct(n).RF is at least an empty struct to fill results
    if ~isfield(DataStruct(n), 'RF'), DataStruct(n).RF = struct(); end

    %% 2.4 Water-level deconvolution
    try
        % (The function signature of makeRFwater_ammon may vary in your codebase)
        [wlr, wlrms, nwl] = makeRFwater_ammon( ...
            R, Z, param.phaseshift, dt, nt, param.waterlevel, param.gauss, param.verbose);

        % Store in DataStruct
        DataStruct(n).RF.wlr    = wlr(:);
        DataStruct(n).RF.wlrms  = wlrms;
        DataStruct(n).RF.nwl    = nwl;
        DataStruct(n).RF.wltime = (dt*(0:nt-1) - param.phaseshift)';

        % Log message
        logMsg = sprintf('[Deconv] Water-level decon OK (gauss=%.2f, wlevel=%.3f)', ...
            param.gauss, param.waterlevel);
        DataStruct = appendHistory(DataStruct, n, logMsg);

    catch ME
        warnMsg = sprintf('[Deconv] Trace %d: Water-level decon failed: %s', n, ME.message);
        warning(warnMsg);
        DataStruct = appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    %% 2.5 Iterative deconvolution
    try
        % (Again, adapt to your actual function signature for iterative decon)
        [itr, itrms] = makeRFitdecon_la_norm( ...
            R, Z, dt, nt, param.phaseshift, param.gauss, param.itmax, param.minderr);
        
        if isempty(itrms)
            % If itrms is empty, no valid solution
            warnMsg = sprintf('[Deconv] Trace %d: empty itrms -> no iterative RF.', n);
            warning(warnMsg);
            DataStruct = appendHistory(DataStruct, n, warnMsg);
            noresult(end+1) = n;
        else
            % Store iterative result
            DataStruct(n).RF.itr    = itr(:);
            DataStruct(n).RF.itrms  = itrms(end);
            DataStruct(n).RF.ittime = (dt*(0:nt-1) - param.phaseshift)';

            logMsg = sprintf('[Deconv] Iterative decon OK (gauss=%.2f, itmax=%d, minderr=%.1e)', ...
                param.gauss, param.itmax, param.minderr);
            DataStruct = appendHistory(DataStruct, n, logMsg);
        end
    catch ME
        warnMsg = sprintf('[Deconv] Trace %d: Iterative decon failed: %s', n, ME.message);
        warning(warnMsg);
        DataStruct = appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end
end

%% 3) Optionally remove traces with no results
if ~isempty(noresult)
    % If desired, remove the invalid or unsuccessful records
    DataStruct(noresult) = [];
    fprintf('Removed %d traces with no valid decon results.\n', length(noresult));
end

elapsedTime = toc(tStart);
disp(['--- Deconvolution completed in ' num2str(elapsedTime,'%.2f') ' s ---']);
end

%% Helper function: appendHistory
function DataStruct = appendHistory(DataStruct, idx, msg)
% APPENDHISTORY  Append a message to DataStruct(idx).ProcHistory.
%                If it does not exist, create it as a cell array.

if ~isfield(DataStruct(idx), 'ProcHistory') || isempty(DataStruct(idx).ProcHistory)
    DataStruct(idx).ProcHistory = {msg};
else
    DataStruct(idx).ProcHistory{end+1} = msg;
end
end