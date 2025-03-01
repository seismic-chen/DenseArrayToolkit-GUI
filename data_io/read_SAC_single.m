function tmpStruct = read_SAC_single(sacFilePath)
% read_SAC  Reads a single Z-component SAC file and assembles a DataStruct array.
%
% Usage:
%   tmpStruct = read_SAC(sacFilePath)
%
% Inputs:
%   sacFilePath : Path of the SAC file
%
% Output:
%   tmpStruct : Structure array, each element contains E/N/Z waveforms
%                and the corresponding station, event, and time-axis info.

% Author: Yunfeng Chen (Refined by ChatGPT)
% Date:   Feb. 2, 2025

%--------------- Read files and fill DataStruct ---------------%
%===  Full path to the Z file and read it ===%
Zfile = fullfile(sacFilePath);
[tZ, dataZ, hdrZ] = fget_sac(Zfile);

%=== Infer the corresponding N/E filenames ===%
Zcomp = deblank(char(hdrZ.stations.kcmpnm)); % channel name from SAC header
Ncomp = strrep(Zcomp,'Z','N');
Ecomp = strrep(Zcomp,'Z','E');
Nfile = strrep(Zfile, Zcomp, Ncomp);
Efile = strrep(Zfile, Zcomp, Ecomp);

%=== Check if N/E files exist; skip if missing ===%
if ~isfile(Nfile) || ~isfile(Efile)
    skipMsg = sprintf('[Warning] Missing N or E file for %s. Skipped.\n', Zfile);
    warning(skipMsg);
    % return an empty struct
    tmpStruct = struct();
    % Optionally record this in ProcHistory
    return;  % Skip this one
end

%=== Read N/E files ===%
[tN, dataN, hdrN] = fget_sac(Nfile);
[tE, dataE, hdrE] = fget_sac(Efile);

tmpStruct = initEmptyStruct();

%------ (A) Time-axis handling ------%
tmpStruct.TimeAxis.dt = hdrZ.times.delta;
tmpStruct.TimeAxis.b  = hdrZ.times.b;
tmpStruct.TimeAxis.e  = hdrZ.times.e;
tmpStruct.TimeAxis.o  = hdrZ.times.o;
if tmpStruct.TimeAxis.o == -12345
    tmpStruct.TimeAxis.o = 0; % SAC often uses -12345 to indicate missing info
end

% Absolute start time (UTC) & time vector
dateNumZ = sacHeaderToDateNum(hdrZ);
tmpStruct.TimeAxis.startTimeUTC = dateNumZ + hdrZ.times.b/86400;
tmpStruct.TimeAxis.endTimeUTC   = dateNumZ + hdrZ.times.e/86400;
tmpStruct.TimeAxis.t            = (0:length(dataZ)-1)' * hdrZ.times.delta + hdrZ.times.b;
tmpStruct.EventInfo.orginTimeUTC= dateNumZ + tmpStruct.TimeAxis.o/86400;

%------ (B) Adjust the lengths of Z/N/E to the smallest length ------%
minLen = min([length(dataZ), length(dataN), length(dataE)]);
if any([length(dataZ), length(dataN), length(dataE)] ~= minLen)
    mismatchMsg = sprintf(['[Info] Length mismatch among Z, N, E. ' ...
        'Truncating to %d.\n'], minLen);
    warning(mismatchMsg);
    tmpStruct.ProcHistory{end+1} = mismatchMsg;
end
dataZ = dataZ(1:minLen);
dataN = dataN(1:minLen);
dataE = dataE(1:minLen);
tmpStruct.TimeAxis.t = tmpStruct.TimeAxis.t(1:minLen);
tmpStruct.TimeAxis.endTimeUTC = tmpStruct.TimeAxis.startTimeUTC + ...
    (minLen - 1)*hdrZ.times.delta/86400;

%------ (C) Store three-component data ------%
tmpStruct.Waveforms.data   = [dataE, dataN, dataZ];
tmpStruct.Waveforms.chName = {Ecomp, Ncomp, Zcomp};

%------ (D) Station and event info ------%
tmpStruct.StationInfo.sta     = deblank(char(hdrZ.station.kstnm));
tmpStruct.StationInfo.stla    = hdrZ.station.stla;
tmpStruct.StationInfo.stlo    = hdrZ.station.stlo;
tmpStruct.StationInfo.stel    = hdrZ.station.stel;
tmpStruct.StationInfo.network = deblank(char(hdrZ.stations.knetwk));

tmpStruct.EventInfo.evla = hdrZ.event.evla;
tmpStruct.EventInfo.evlo = hdrZ.event.evlo;
tmpStruct.EventInfo.evdp = hdrZ.event.evdp;
tmpStruct.EventInfo.mag  = hdrZ.event.mag;
tmpStruct.EventInfo.evid = datestr(tmpStruct.EventInfo.orginTimeUTC,'yyyymmddHHMMSS');

%------ (E) File header info and initial processing history ------%
tmpStruct.Header.filenameZ = Zfile;
tmpStruct.Header.filenameN = Nfile;
tmpStruct.Header.filenameE = Efile;
initLog = sprintf('Data read from SAC (E, N, Z) on %s', datestr(now, 31));
tmpStruct.ProcHistory{end+1} = initLog;

%------ (F) Check for inconsistencies among Z/N/E headers (example) ------%
[inconsistency, detailMsg] = checkHeaderInconsistency(hdrZ, hdrN, hdrE);
if inconsistency
    warning(detailMsg);
    tmpStruct.ProcHistory{end+1} = detailMsg;
end

end

%% ---------------------- Helper functions below ---------------------- %%
function tmpStruct = initEmptyStruct()
% initEmptyStruct  Creates an empty structure template for pre-allocation.
tmpStruct = struct( ...
    'Waveforms',   [], ...  % waveform data: Nt x Nch
    'TimeAxis',    [], ...  % absolute time, sampling rate, etc.
    'StationInfo', [], ...  % station metadata
    'EventInfo',   [], ...  % event metadata
    'Header',      [], ...  % file path / header info
    'RF',          [], ...  % receiver function (optional)
    'TravelInfo',  [], ...  % travel time / azimuth info (optional)
    'ProcHistory', {{}} ... % cell array for processing history
    );
end

function dateNum = sacHeaderToDateNum(hdr)
% sacHeaderToDateNum  Converts SAC header's yday-based time fields to a MATLAB datenum.
%
% It requires the following fields in hdr:
%   hdr.event.nzyear, hdr.event.nzjday, hdr.event.nzhour,
%   hdr.event.nzmin,  hdr.event.nzsec,  hdr.event.nzmsec
%
% If any are -12345 (often meaning missing data in SAC), an error is raised.

nzyear = hdr.event.nzyear;
nzjday = hdr.event.nzjday;
nzhour = hdr.event.nzhour;
nzmin  = hdr.event.nzmin;
nzsec  = hdr.event.nzsec;
nzmsec = hdr.event.nzmsec;

if any([nzyear, nzjday, nzhour, nzmin, nzsec] == -12345)
    error('Invalid SAC header: missing time info (year/day/hour/min/sec).');
end

[nzday, nzmonth] = jday_to_day(nzyear, nzjday);
date_str = sprintf('%04d-%02d-%02d %02d:%02d:%02d', ...
    nzyear, nzmonth, nzday, nzhour, nzmin, nzsec);
dateNum = datenum(date_str, 'yyyy-mm-dd HH:MM:SS');
dateNum = dateNum + (nzmsec/1000)/86400; % Add milliseconds
end

function [inconsistency, detailMsg] = checkHeaderInconsistency(hdrZ, hdrN, hdrE)
% checkHeaderInconsistency  Example function to detect header differences across Z/N/E.
% You can compare station coords, sample rates, etc. and report any major discrepancies.

inconsistency = false;
detailMsg = '';

% Example: compare station lat among Z, N, E
if abs(hdrZ.station.stla - hdrN.station.stla) > 1e-5 || ...
        abs(hdrZ.station.stla - hdrE.station.stla) > 1e-5
    inconsistency = true;
    detailMsg = sprintf('[Warning] Station lat mismatch: Z=%.6f, N=%.6f, E=%.6f\n',...
        hdrZ.station.stla, hdrN.station.stla, hdrE.station.stla);
end

% You can add more comparisons as needed (e.g., station.stlo, sampling rates, etc.)
end