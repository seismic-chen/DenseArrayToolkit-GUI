function [CommonEventGather,matchIndex] = getCommonEventGather(DataStruct, EventID)
% GETCOMMANDEVENTGATHER 提取来自同一事件的所有台站数据
%
% 使用方法:
%   CommonEventGather = getCommonEventGather(DataStruct, EventID)
%
% 输入参数:
%   DataStruct - struct 数组，包含字段:
%       .EventInfo.evid - 事件标识符（字符串）
%   EventID    - 要提取的事件标识符（字符串）
%
% 输出参数:
%   CommonEventGather - struct 数组，仅包含与 EventID 匹配的记录
%
% 示例:
%   CommonGather = getCommonEventGather(DataStruct, 'EV12345');

    %% 1. 输入参数验证
    if nargin < 2
        error('getCommonEventGather:InsufficientInputs', '需要两个输入参数: DataStruct 和 EventID.');
    end
    
    if ~isstruct(DataStruct)
        error('getCommonEventGather:InvalidDataStruct', 'DataStruct 必须是一个结构数组.');
    end
    
    if ~isfield(DataStruct, 'EventInfo') || ~isfield(DataStruct(1).EventInfo, 'evid')
        error('getCommonEventGather:MissingField', 'DataStruct 中必须包含 EventInfo.evid 字段.');
    end
    
    if ~(ischar(EventID) || isstring(EventID))
        error('getCommonEventGather:InvalidEventID', 'EventID 必须是字符串或字符数组.');
    end

    %% 2. 使用逻辑索引筛选匹配的记录
    % 使用 arrayfun 创建逻辑数组，检查每个记录的 EventID 是否匹配
    isMatch = arrayfun(@(x) strcmp(x.EventInfo.evid, EventID), DataStruct);
    matchIndex = find(isMatch);
    % 根据逻辑数组提取匹配的记录
    CommonEventGather = DataStruct(isMatch);
    %% 3. 处理未找到的情况
    if isempty(CommonEventGather)
        warning('getCommonEventGather:NoMatch', '未找到与 EventID "%s" 匹配的记录.', EventID);
    else
        fprintf('找到 %d 条与 EventID "%s" 匹配的记录.\n', length(CommonEventGather), EventID);
    end
end