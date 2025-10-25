function [CommonStationGather, matchIndex] = getCommonStationGather(DataStruct, station, matchOptions)
% GETCOMMONSTATIONGATHER Extract all records from specified station(s)
%
% Syntax:
%   [CommonStationGather, matchIndex] = getCommonStationGather(DataStruct, station, matchOptions)
%
% Input Parameters:
%   DataStruct   - struct array containing fields:
%                   .StationInfo.sta - Station name (string)
%   station      - Station name(s) to extract (string, char, or string array)
%   matchOptions - (optional) Matching options structure with fields:
%                   .CaseInsensitive - Case-insensitive matching (logical, default: false)
%                   .PartialMatch    - Partial string matching (logical, default: false)
%
% Output Parameters:
%   CommonStationGather - struct array containing only records matching specified station(s)
%   matchIndex          - indices of matching records in original DataStruct
%
% Examples:
%   % Extract single station, exact match, case-sensitive
%   [Gather, idx] = getCommonStationGather(DataStruct, 'STA1');
%
%   % Extract multiple stations, case-insensitive, exact match
%   [Gather, idx] = getCommonStationGather(DataStruct, {'sta1', 'Sta2'}, struct('CaseInsensitive', true));
%
%   % Extract stations containing 'STA', case-insensitive, partial match
%   [Gather, idx] = getCommonStationGather(DataStruct, 'STA', struct('CaseInsensitive', true, 'PartialMatch', true));
%

    %% 1. Input Validation
    if nargin < 2
        error('getCommonStationGather:InsufficientInputs', 'At least two input arguments required: DataStruct and station.');
    end

    if ~isstruct(DataStruct)
        error('getCommonStationGather:InvalidDataStruct', 'DataStruct must be a structure');
    end

    % Check if DataStruct contains StationInfo and StationInfo.sta
    if ~all(isfield(DataStruct, 'StationInfo')) || ~all(isfield([DataStruct.StationInfo], 'sta'))
        error('getCommonStationGather:MissingField', 'Each element in DataStruct must contain StationInfo.sta field.');
    end

    % Check if station is string or string array
    if ~(ischar(station) || isstring(station) || (iscellstr(station) || isstring(station)))
        error('getCommonStationGather:InvalidStation', 'station must be a string, char array, or cell array of strings.');
    end

    % Handle matchOptions
    if nargin < 3 || isempty(matchOptions)
        matchOptions = struct();
    end

    % Set default options
    if ~isfield(matchOptions, 'CaseInsensitive')
        matchOptions.CaseInsensitive = false;
    end
    if ~isfield(matchOptions, 'PartialMatch')
        matchOptions.PartialMatch = false;
    end

    %% 2. Extract all station names
    stationNames = {};
    for n = 1:length(DataStruct)
        stationNames{end+1} = DataStruct(n).StationInfo.sta;
    end
    %% 3. Generate matching logical array
    if ischar(station) || isstring(station)
        station = {station};  % Convert to cell array
    elseif iscell(station)
        if ~all(cellfun(@(x) ischar(x) || isstring(x), station))
            error('getCommonStationGather:InvalidStation', 'All station cell elements must be strings or char');
        end
    end

    % Initialize matching logical array
    isMatch = false(size(stationNames));

    for i = 1:length(station)
        currentStation = station{i};
        if matchOptions.CaseInsensitive
            currentStation = lower(currentStation);
            stationNamesLower = lower(stationNames);
        else
            stationNamesLower = stationNames;
        end

        if matchOptions.PartialMatch
            isMatch = isMatch | contains(stationNamesLower, currentStation);
        else
            isMatch = isMatch | strcmp(stationNamesLower, currentStation);
        end
    end

    %% 4. Extract matching records and indices
    matchIndex = find(isMatch);
    CommonStationGather = DataStruct(isMatch);

    %% 5. Handle no matches found
    if isempty(CommonStationGather)
         warning('getCommonStationGather:NoMatch', 'No records found matching the specified station(s).');
    else
        fprintf('Found %d records matching the specified station \n', length(CommonStationGather));
    end
end