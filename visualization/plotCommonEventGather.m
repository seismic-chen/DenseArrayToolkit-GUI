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

    %% 1) check input
    if nargin < 3
        warning('plotSingleEventRFs requires at least 3 inputs: DataStruct, eventID, plot_type. Use trace index as default.');
        axis_type = 'trace';
        plot_type = 'imagesc';
    end

    if nargin < 2 || isempty(EventID)
        EventID = 'all';
    end

    stmn = {};
    for i = 1:numel(DataStruct)
        stmn{i} = [DataStruct(i).StationInfo.sta];
    end
    nsta = length(sort(unique(stmn)));
    stnet = strsplit(stmn{1});
    if strcmp(stnet{1}(1:2),'SL')==1
        nsta = nsta+3;
    end

    %% 2) 
    rfmatrix = zeros(length(DataStruct(1).RF.itr),nsta);
    distArr  = [];
    foundIndices = [];

    if strcmp(EventID, 'all')
        % 
        foundIndices = 1:numel(DataStruct);
    else
        % 
        isMatch = arrayfun(@(x) strcmp(x.EventInfo.evid, EventID), DataStruct);
        foundIndices = find(isMatch);
    end

    if isempty(foundIndices)
        warning('No data found for eventID = %s', EventID);
        return;
    end

    for idx = foundIndices
        if isfield(DataStruct(idx), 'RF') && isfield(DataStruct(idx).RF,'itr') ...
                && ~isempty(DataStruct(idx).RF.itr)
                      % [Nt x 1] -> add as new column
            stanm = DataStruct(idx).StationInfo.sta;
            tmp = str2double(strsplit(stanm,stanm(1:2)));
            stanum = tmp(2);
            rfmatrix(:, stanum) = DataStruct(idx).RF.itr;

        else
            continue;
        end

        if isfield(DataStruct(idx), 'TravelInfo') && isfield(DataStruct(idx).TravelInfo, 'distDeg')
            distArr(end+1) = DataStruct(idx).TravelInfo.distDeg; % distance in deg
        else
            distArr(end+1) = NaN;  
        end
    end

    if isempty(rfmatrix)
        warning('Found %d records, but none has valid RF.itr data.', ...
                 numel(foundIndices));
        return;
    end

    if strcmp(stnet{1}(1:2),'SL')==1
        rfmatrix = fliplr(rfmatrix);
    end


    %% 3) time axis
    %   
    firstIdx = foundIndices(1);
    if isfield(DataStruct(firstIdx), 'RF') && isfield(DataStruct(firstIdx).RF,'ittime') ...
            && ~isempty(DataStruct(firstIdx).RF.ittime)
        t = DataStruct(firstIdx).RF.ittime;
    else
        warning('No valid time axis (RF.ittime) found in the first matched record.');
        return;
    end
   
    %% 4) plot
    if nargin < 5 || isempty(UIAxes)
        % 
        figure('Name',sprintf('RF for %s',EventID),'Color','white',...
               'Position',[10 10 900 350]);
        ax = gca;
    else
        % 
        ax = UIAxes;
        axes(ax); 
        cla(ax);  
        hold(ax, 'on'); 
    end

    switch axis_type
        case 'distance'
            xrange = distArr;
            xlabelstr = 'Distance (deg)';
        case 'trace'
            xrange = 1:nsta;
            xlabelstr = '# of station';  
    end
    switch plot_type
        case 'wigb'
            wigb(rfmatrix, 2, xrange, t);
            xlabel(ax, xlabelstr);
            xlim(ax,[min(xrange) max(xrange)])
        case 'imagesc'
            imagesc(xrange,t,rfmatrix,'Parent', ax);
            % colormap(ax,seismic(1))
            colormap(ax,'gray')
            clim(ax,[-0.1 0.1])
            xlabel(ax, xlabelstr);
            xlim(ax,[min(xrange) max(xrange)])
    end
    ylim(ax, [t(1), 20]);    
    set(ax,'YDir','reverse')
    ylabel(ax, 'Time (sec)');
    title(ax,['RF (ID: ',EventID,')'])
    set(ax, 'FontSize',18, 'LineWidth',1.5, 'XMinorTick','on');

    if strcmp(EventID, 'all')
        title(ax, sprintf('All Events : showing %d station(s)', size(rfmatrix,2)));
    else
        % title(ax, sprintf('Event %s : showing %d station(s)', EventID, size(rfmatrix,2)));
    end
end