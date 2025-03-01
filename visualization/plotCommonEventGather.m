function plotCommonEventGather(DataStruct, EventID, axis_type, plot_type, UIAxes)
% PLOTSINGLEEVENTRFS - Plot receiver functions for all stations of a single event or all events
%
% Usage:
%   plotSingleEventRFs(DataStruct, eventID, plot_type, UIAxes)
%
% Inputs:
%   DataStruct : struct array, each element has fields:
%       .EventInfo.evid           : event ID (string)
%       .RF.ittime (time axis)    : e.g. [Nt x 1]
%       .RF.itr (receiver func)   : e.g. [Nt x 1]
%       .TravelInfo.distDeg       : distance in degrees
%   eventID    : (optional) the event ID you want to plot. If not provided, plot all events.
%   axis_type  : 'distance' (use distDeg as X-axis) or 'trace' (use trace index)
%   plot_type  : 'wigb' to use wiggle plot or 'imagesc' to use image plot
%   UIAxes     : (optional) Axes handle to plot in. If not provided, a new figure is created.
%
% Example:
%   plotSingleEventRFs(DataStruct, 'EV12345678', 'distance', UIAxes);
%   plotSingleEventRFs(DataStruct, [], 'trace'); % Plot all events

    %% 1) 输入检查
    if nargin < 3
        warning('plotSingleEventRFs requires at least 3 inputs: DataStruct, eventID, plot_type. Use trace index as default.');
        axis_type = 'trace';
        plot_type = 'imagesc';
    end

    if nargin < 2 || isempty(EventID)
        % 如果没有提供 EventID，则绘制所有波形
        EventID = 'all';
    end

    %% 2) 收集满足条件的 RF 波形
    rfmatrix = [];
    distArr  = [];
    foundIndices = [];

    if strcmp(EventID, 'all')
        % 如果 EventID 是 'all'，则收集所有波形
        foundIndices = 1:numel(DataStruct);
    else
        % 否则，收集匹配指定 EventID 的波形
        isMatch = arrayfun(@(x) strcmp(x.EventInfo.evid, EventID), DataStruct);
        foundIndices = find(isMatch);
    end

    if isempty(foundIndices)
        warning('No data found for eventID = %s', EventID);
        return;
    end

    for idx = foundIndices
        % 确保 RF.itr 存在且非空
        if isfield(DataStruct(idx), 'RF') && isfield(DataStruct(idx).RF,'itr') ...
                && ~isempty(DataStruct(idx).RF.itr)
            rfmatrix(:, end+1) = DataStruct(idx).RF.itr;          % [Nt x 1] -> add as new column
        else
            % 若这个记录没有 RF.itr, 跳过
            continue;
        end

        if isfield(DataStruct(idx), 'TravelInfo') && isfield(DataStruct(idx).TravelInfo, 'distDeg')
            distArr(end+1) = DataStruct(idx).TravelInfo.distDeg; % distance in deg
        else
            distArr(end+1) = NaN;  % 或者赋予一个默认值
        end
    end

    % 若最终 rfmatrix 为空, 说明虽然找到匹配事件, 但没有有效RF波形
    if isempty(rfmatrix)
        warning('Found %d records, but none has valid RF.itr data.', ...
                 numel(foundIndices));
        return;
    end

    %% 3) 确定时间轴
    %   从第一个有 RF.ittime 的记录中获取
    firstIdx = foundIndices(1);
    if isfield(DataStruct(firstIdx), 'RF') && isfield(DataStruct(firstIdx).RF,'ittime') ...
            && ~isempty(DataStruct(firstIdx).RF.ittime)
        t = DataStruct(firstIdx).RF.ittime;
    else
        warning('No valid time axis (RF.ittime) found in the first matched record.');
        return;
    end

    [nt, nx] = size(rfmatrix);
    
    %% 4) 作图
    if nargin < 5 || isempty(UIAxes)
        % 如果没有提供 UIAxes，则新建一个 figure
        figure('Name',sprintf('RF for %s',EventID),'Color','white',...
               'Position',[200 200 1200 700]);
        ax = gca;
    else
        % 如果提供了 UIAxes，则在该 Axes 中画图
        ax = UIAxes;
        axes(ax); % 显式设置当前 Axes
        cla(ax);  % 清空现有内容
        hold(ax, 'on'); % 确保后续绘图在同一 Axes
    end

    switch axis_type
        case 'distance'
            xrange = distArr;
            xlabelstr = 'Distance (deg)';
        case 'trace'
            xrange = 1:nx;
            xlabelstr = 'Trace index';  
    end
    switch plot_type
        case 'wigb'
            wigb(rfmatrix, 2, xrange, t);
            xlabel(ax, xlabelstr);
            xlim(ax,[min(xrange) max(xrange)])
        case 'imagesc'
            imagesc(xrange,t,rfmatrix,'Parent', ax);
            colormap(ax,seismic(1))
            caxis(ax,[-0.1 0.1])
            xlabel(ax, xlabelstr);
            xlim(ax,[min(xrange) max(xrange)])
    end
    ylim(ax, [t(1), 20]);    % 只显示 0~20s 区间，可按需调整
    set(ax,'YDir','reverse')
    ylabel(ax, 'Time (sec)');
    set(ax, 'FontSize',14, 'LineWidth',1, 'XMinorTick','on');

    if strcmp(EventID, 'all')
        title(ax, sprintf('All Events : showing %d station(s)', size(rfmatrix,2)));
    else
        title(ax, sprintf('Event %s : showing %d station(s)', EventID, size(rfmatrix,2)));
    end
end