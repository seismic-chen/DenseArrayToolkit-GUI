function [CommonEventGather,matchIndex] = getCommonEventGather(DataStruct, EventID)
% GETCOMMONEVENTGATHER Extract station data from the same seismic event
%
% Syntax:
%   CommonEventGather = getCommonEventGather(DataStruct, EventID)
%
% Input Parameters:
%   DataStruct - struct array containing fields:
%       .EventInfo.evid - Event identifier (string)
%   EventID    - Target event identifier to extract (string)
%
% Output Parameters:
%   CommonEventGather - struct array containing only records matching EventID
%   matchIndex        - indices of matching records in the original DataStruct
%
% Example:
%   CommonGather = getCommonEventGather(DataStruct, 'EV12345');

    %% 1. 输入参数验证
    if nargin < 2
        error('getCommonEventGather:InsufficientInputs',...
            'Two input arguments required: DataStruct and EventID');
    end
    
    if ~isstruct(DataStruct)
        error('getCommonEventGather:InvalidDataStruct', 'DataStruct must be a structure variable');
    end
    
    if ~isfield(DataStruct, 'EventInfo') || ~isfield(DataStruct(1).EventInfo, 'evid')
        error('getCommonEventGather:MissingField', 'must contain EventInfo.evid field');
    end
    
    if ~(ischar(EventID) || isstring(EventID))
        error('getCommonEventGather:InvalidEventID', 'EventID must be a string or character');
    end

    %% 2. Filter matching records using logical indexing
    % Use arrayfun to create logical array checking EventID match for each record
    isMatch = arrayfun(@(x) strcmp(x.EventInfo.evid, EventID), DataStruct);
    matchIndex = find(isMatch);
    % Extract matching records based on logical array
    CommonEventGather = DataStruct(isMatch);
    %% 3. Handle no matches found
    if isempty(CommonEventGather)
        warning('getCommonEventGather:NoMatch', 'No records found matching EventID "%s".', EventID);
    else
        fprintf('Found %d records matching EventID: "%s"\n', length(CommonEventGather), EventID);
    end
end