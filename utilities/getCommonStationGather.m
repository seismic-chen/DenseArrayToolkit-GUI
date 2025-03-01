function [CommonStationGather, matchIndex] = getCommonStationGather(DataStruct, station, matchOptions)
% GETCOMMONSTATIONGATHER 提取来自指定台站的所有记录
%
% 使用方法:
%   [CommonStationGather, matchIndex] = getCommonStationGather(DataStruct, station, matchOptions)
%
% 输入参数:
%   DataStruct   - struct 数组，包含字段:
%                   .StationInfo.sta - 台站名称（字符串）
%   station      - 要提取的台站名称（字符串）或台站名称的字符串数组
%   matchOptions - (可选) 匹配选项，结构体包含字段:
%                   .CaseInsensitive - 是否忽略大小写（逻辑值，默认: false）
%                   .PartialMatch     - 是否进行部分匹配（逻辑值，默认: false）
%
% 输出参数:
%   CommonStationGather - struct 数组，仅包含与指定台站匹配的记录
%   matchIndex          - 数组，匹配记录在原 DataStruct 中的索引
%
% 示例:
%   % 提取单个台站，精确匹配，区分大小写
%   [Gather, idx] = getCommonStationGather(DataStruct, 'STA1');
%
%   % 提取多个台站，忽略大小写，精确匹配
%   [Gather, idx] = getCommonStationGather(DataStruct, {'sta1', 'Sta2'}, struct('CaseInsensitive', true));
%
%   % 提取包含 'STA' 的台站，忽略大小写，部分匹配
%   [Gather, idx] = getCommonStationGather(DataStruct, 'STA', struct('CaseInsensitive', true, 'PartialMatch', true));
%
% 作者: <你的名字>
% 日期: <日期>

    %% 1. 输入参数验证
    if nargin < 2
        error('getCommonStationGather:InsufficientInputs', '需要至少两个输入参数: DataStruct 和 station.');
    end

    if ~isstruct(DataStruct)
        error('getCommonStationGather:InvalidDataStruct', 'DataStruct 必须是一个结构数组.');
    end

    % 检查 DataStruct 是否包含 StationInfo 和 StationInfo.sta
    if ~all(isfield(DataStruct, 'StationInfo')) || ~all(isfield([DataStruct.StationInfo], 'sta'))
        error('getCommonStationGather:MissingField', 'DataStruct 中的每个元素必须包含 StationInfo.sta 字段.');
    end

    % 检查 station 是否为字符串或字符串数组
    if ~(ischar(station) || isstring(station) || (iscellstr(station) || isstring(station)))
        error('getCommonStationGather:InvalidStation', 'station 必须是字符串、字符串数组或字符数组的单元格.');
    end

    % 处理 matchOptions
    if nargin < 3 || isempty(matchOptions)
        matchOptions = struct();
    end

    % 设置默认匹配选项
    if ~isfield(matchOptions, 'CaseInsensitive')
        matchOptions.CaseInsensitive = false;
    end
    if ~isfield(matchOptions, 'PartialMatch')
        matchOptions.PartialMatch = false;
    end

    %% 2. 提取所有台站名称
    stationNames = {};
    for n = 1:length(DataStruct)
        stationNames{end+1} = DataStruct(n).StationInfo.sta;
    end
    %% 3. 生成匹配逻辑数组
    if ischar(station) || isstring(station)
        station = {station};  % 转换为单元格数组
    elseif iscell(station)
        % 确保所有元素都是字符串
        if ~all(cellfun(@(x) ischar(x) || isstring(x), station))
            error('getCommonStationGather:InvalidStation', 'station 单元格数组中的所有元素必须是字符串或字符数组.');
        end
    end

    % 初始化匹配逻辑数组
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

    %% 4. 提取匹配的记录和索引
    matchIndex = find(isMatch);
    CommonStationGather = DataStruct(isMatch);

    %% 5. 处理未找到的情况
    if isempty(CommonStationGather)
        warning('getCommonStationGather:NoMatch', '未找到与指定台站匹配的记录.');
    else
        fprintf('找到 %d 条与指定台站匹配的记录.\n', length(CommonStationGather));
    end
end