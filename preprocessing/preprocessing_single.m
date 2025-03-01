function [DataStruct,removeFlag] = preprocessing_single(DataStruct, param)
% PREPROCESSING  Applies a set of preprocessing steps to DataStruct,
%                including:
%                - Distance filtering
%                - P arrival calculation
%                - SNR computation
%                - DC/trend removal
%                - ENZ -> TRZ rotation
%                - Bandpass filtering
%                - Windowing around the P wave
%                - Resampling
%
% Usage:
%   DataStruct = preprocessing(DataStruct, param)
%
% Inputs:
%   DataStruct : Struct array, each element typically corresponds to one
%                three-component record (E, N, Z). Required fields include:
%       - Waveforms.data      (raw data, [Nt x 3])
%       - TimeAxis.b, TimeAxis.o, TimeAxis.dt (from SAC header)
%       - EventInfo.evla, evlo, evdp
%       - StationInfo.stla, stlo
%   param : Struct of preprocessing parameters. Default fields:
%       param.tstart         (start time offset, not always used)
%       param.tend           (end time offset, not always used)
%       param.sig_leader     (time to pad before P wave)
%       param.record_len     (length of desired window after sig_leader)
%       param.lows, highs    (bandpass filter corners, e.g. [0.1, 1.2] Hz)
%       param.resample_period (new sample period, e.g. 0.1s => 10 Hz)
%
% Outputs:
%   DataStruct : Updated struct array with fields:
%       - Waveforms.dataProcessed
%       - Waveforms.chNameProcessed
%       - TimeAxis.t_resample
%       - TimeAxis.dt_resample
%       - RF.snr  (signal-to-noise ratio on Z)
%       - TravelInfo (distance, azimuth, P arrival time, etc.)
%       - ProcHistory (new logs appended)
%
% Example:
%   config.PreprocessingParam.lows = 0.1;
%   config.PreprocessingParam.highs = 1.2;
%   ...
%   DataStruct = preprocessing(DataStruct, config.PreprocessingParam);
%
% Author: Yunfeng Chen
% Date  : Jan 9, 2025

%% 1. Parameter defaults / initialization
if nargin < 2, param = struct(); end

% Check for missing fields, set defaults:
if ~isfield(param, 'tstart'),          param.tstart = -60; end
if ~isfield(param, 'tend'),            param.tend = 600;   end
if ~isfield(param, 'sig_leader'),      param.sig_leader = 30; end
if ~isfield(param, 'record_len'),      param.record_len = 120; end
if ~isfield(param, 'lows'),            param.lows = 0.1;  end
if ~isfield(param, 'highs'),           param.highs = 1.2; end
if ~isfield(param, 'resample_period'), param.resample_period = 0.1; end

% This vector will collect indices of records that fail certain checks
removeFlag = false;
%% 2. Main processing

%% 2.1 Basic checks: is waveform data present?
if ~isfield(DataStruct.Waveforms, 'data') || isempty(DataStruct.Waveforms.data)
    msgSkip = '[Warning] No waveform data found -> skipping this record.';
    DataStruct.ProcHistory{end+1} = msgSkip;
    removeFlag = true;
    return;
end

seis = DataStruct.Waveforms.data;  % [Nt x 3], typically ENZ
[nt, nchan] = size(seis);

% Check if we have at least 3 components
if nchan < 3
    msgSkip = '[Warning] Fewer than 3 components -> skipping.';
    DataStruct.ProcHistory{end+1} = msgSkip;
    removeFlag = true;
    return;
end

%% 2.2 Retrieve event and station info
evla = DataStruct.EventInfo.evla;
evlo = DataStruct.EventInfo.evlo;
evdp = DataStruct.EventInfo.evdp;
stla = DataStruct.StationInfo.stla;
stlo = DataStruct.StationInfo.stlo;

% Compute distance, azimuth using haversine formula
try
    [az, baz, degree] = distHaversine(evla, evlo, stla, stlo);
catch ME
    msgErr = sprintf('[Error] distHaversine failed: %s', ME.message);
    DataStruct.ProcHistory{end+1} = msgErr;
    removeFlag = true;
    return;
end

% Store results into TravelInfo
DataStruct.TravelInfo.distDeg = degree;
DataStruct.TravelInfo.az      = az;
DataStruct.TravelInfo.baz     = baz;

% Keep data only if epicentral distance is in [30, 95]
if degree < 30 || degree > 95
    msgRange = sprintf('[Info] dist=%.1f out of [30,95] -> removed.', degree);
    DataStruct.ProcHistory{end+1} = msgRange;
    removeFlag = true;
    return;
end

%% 2.3 Compute P-wave arrival time via taupTime
try
    ptime = taupTime('prem', evdp, 'P', 'deg', degree);
    if isempty(ptime)
        msgPT = '[Warning] No P-wave arrival found -> removing.';
        DataStruct.ProcHistory{end+1} = msgPT;
        removeFlag = true;
        return;
    end
    pTime = ptime(1).time;

    % Save to DataStruct
    DataStruct.TravelInfo.pTime    = pTime;
    DataStruct.TravelInfo.rayParam = ptime(1).rayParam;
catch ME
    msgPT = sprintf('[Error] taupTime: %s', ME.message);
    DataStruct.ProcHistory{end+1} = msgPT;
    removeFlag = true;
    return;
end

%% 2.4 Compute SNR (using a separate function calcSNR)
% Attempt to retrieve dt, b, o from TimeAxis
if isfield(DataStruct.TimeAxis, 'dt'), dt = DataStruct.TimeAxis.dt; else, dt = []; end
if isfield(DataStruct.TimeAxis, 'b'),  b  = DataStruct.TimeAxis.b;  else, b  = 0;  end
if isfield(DataStruct.TimeAxis, 'o'),  o  = DataStruct.TimeAxis.o;  else, o  = 0;  end
if isfield(DataStruct.TimeAxis, 't'),  t  = DataStruct.TimeAxis.t;  else, t  = (0:nt-1)*dt; end

Z = seis(:,3);  % Z is typically the 3rd column
try
    snrVal = calcSNR(Z, dt, b, o, pTime);
    DataStruct.RF.snr = snrVal;
    DataStruct.ProcHistory{end+1} = ...
        sprintf('[Info] SNR computed on Z: %.2f', snrVal);
catch ME
    msgSNR = sprintf('[Warning] calcSNR failed: %s', ME.message);
    DataStruct.ProcHistory{end+1} = msgSNR;
    DataStruct.RF.snr = -999;  % some flag indicating invalid SNR
end

%% 2.5 Remove DC offset & linear trend
seis = removeSeisDC(seis);
seis = removeSeisTrend(seis);
DataStruct.ProcHistory{end+1} = '[Info] DC & trend removed.';

%% 2.6 Rotate from ENZ to TRZ
try
    [seis_rot, ~] = rotateSeisENZtoTRZ(seis, baz);
    DataStruct.ProcHistory{end+1} = ...
        sprintf('[Info] rotateSeisENZtoTRZ done, baz=%.2f.', baz);
catch ME
    msgRot = sprintf('[Error] rotateSeisENZtoTRZ failed: %s', ME.message);
    DataStruct.ProcHistory{end+1} = msgRot;
    removeFlag = true;
    return;
end

%% 2.7 Taper and filter
seis_taper = taperSeis(seis_rot, 0.2);
DataStruct.ProcHistory{end+1} = '[Info] taperSeis applied.';

if ~isempty(dt)
    % Apply a bandpass filter (e.g. 3rd-order Butterworth)
    seis_flt = bandpassSeis(seis_taper, dt, param.lows, param.highs, 3);
    DataStruct.ProcHistory{end+1} = ...
        sprintf('[Info] bandpassSeis (%.2f-%.2f Hz) done.', param.lows, param.highs);
else
    seis_flt = seis_taper;
    DataStruct.ProcHistory{end+1} = ...
        '[Warning] dt not found -> skipping bandpassSeis.';
end

%% 2.8 Chop the data window around the P wave
%   param.sig_leader = time before P (e.g. 30s)
%   param.record_len = total window length after sig_leader
startCut  = pTime - param.sig_leader;
finishCut = startCut + param.sig_leader + param.record_len;

[seis_cut, t_cut] = chopSeis(seis_flt, t, startCut, finishCut);
DataStruct.ProcHistory{end+1} = ...
    sprintf('[Info] chopSeis around P=%.1f sec: [%.1f -> %.1f].', ...
    pTime, startCut, finishCut);

%% 2.9 Resample to a new sampling rate
[seis_rsp, dt_rsp, t_rsp] = resampleSeis(seis_cut, t_cut, param.resample_period);
DataStruct.ProcHistory{end+1} = ...
    sprintf('[Info] resampled from dt=%.3f s to %.3f s.', dt, dt_rsp);

%% 2.10 Save the processed data in DataStruct
DataStruct.Waveforms.dataProcessed   = seis_rsp;
DataStruct.Waveforms.chNameProcessed = {'T','R','Z'};  % after rotation
DataStruct.TimeAxis.t_resample       = t_rsp;
DataStruct.TimeAxis.dt_resample      = dt_rsp;
end

%% ===================================================================== %%
function snrVal = calcSNR(Z, dt, b, o, p_time)
% calcSNR  Computes a simple SNR for the Z component as a ratio of
%          variance(signal) / variance(noise).
%
% Usage:
%   snrVal = calcSNR(Z, dt, b, o, p_time)
%
% Inputs:
%   Z      : [Nt x 1] vertical waveform
%   dt     : sampling interval (seconds)
%   b, o   : from the SAC header (b = start time, o = origin time)
%   p_time : P-wave arrival time (seconds, relative to origin 'o')
%
% The noise window is defined as [p_time - 105, p_time - 5]
% The signal window is [p_time - 5, p_time + 5]
%
% Output:
%   snrVal : ratio of var(signal) / var(noise). Larger = better.
%
% Author:  Yunfeng Chen
% Date  :  Feb. 3, 2025

if isempty(dt) || dt <= 0
    error('calcSNR:InvalidDT','dt must be a positive number.');
end

% Define the time windows
noiseStart = p_time - 105;
noiseEnd   = p_time -   5;
sigStart   = p_time -   5;
sigEnd     = p_time +   5;

N = length(Z);

% Convert times to sample indices:
%   idx = round( (time + (o - b)) / dt );
% This accounts for the difference between the record start time (b)
% and the event origin time (o).
idxNoise1 = round((noiseStart + (o - b)) / dt);
idxNoise2 = round((noiseEnd   + (o - b)) / dt);
idxSig1   = round((sigStart  + (o - b)) / dt);
idxSig2   = round((sigEnd    + (o - b)) / dt);

% Bound checks to ensure valid indices
idxNoise1 = max(idxNoise1, 1);
idxNoise2 = min(idxNoise2, N);
idxSig1   = max(idxSig1,   1);
idxSig2   = min(idxSig2,   N);

% If the noise or signal windows are invalid or empty, raise an error
if idxNoise1 >= idxNoise2 || idxSig1 >= idxSig2
    error('calcSNR:InvalidWindow', ...
        'Noise or signal window indices are invalid or overlapping.');
end

% Extract noise and signal segments
noi = Z(idxNoise1 : idxNoise2);
sig = Z(idxSig1   : idxSig2);

% SNR: ratio of variance
snrVal = var(sig) / var(noi);
end