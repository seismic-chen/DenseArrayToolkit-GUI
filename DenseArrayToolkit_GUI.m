classdef DenseArrayToolkit_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        FileMenu                      matlab.ui.container.Menu
        LoadMenu                      matlab.ui.container.Menu
        SaveMenu                      matlab.ui.container.Menu
        ToolsMenu                     matlab.ui.container.Menu
        RadonTransformMenu            matlab.ui.container.Menu
        RankReductionMenu             matlab.ui.container.Menu
        ImagingMenu                   matlab.ui.container.Menu
        HkStackingMenu                matlab.ui.container.Menu
        CCPMenu                       matlab.ui.container.Menu
        Migration2DMenu               matlab.ui.container.Menu
        Migration3DMenu               matlab.ui.container.Menu
        EditingMenu                   matlab.ui.container.Menu
        ViewMenu                      matlab.ui.container.Menu
        InfoTabGroup                  matlab.ui.container.TabGroup
        EventTab                      matlab.ui.container.Tab
        EventTable                    matlab.ui.control.Table
        StationTab                    matlab.ui.container.Tab
        StationTable                  matlab.ui.control.Table
        ButtonGroup                   matlab.ui.container.ButtonGroup
        StationMapPanel               matlab.ui.container.Panel
        EventMapPanel                 matlab.ui.container.Panel
        statusLabel                   matlab.ui.control.Label
        processingWorkflowLabel       matlab.ui.control.Label
        setupPathsButton              matlab.ui.control.Button
        loadConfigButton              matlab.ui.control.Button
        plotStationsEventsButton      matlab.ui.control.Button
        plotCommonEventGatherButton   matlab.ui.control.Button
        plotComonStationGatherButton  matlab.ui.control.Button
        visualizationOptionLabel      matlab.ui.control.Label
        plotWaveformsButton           matlab.ui.control.Button
        stackButton                   matlab.ui.control.Button
        deconvButton                  matlab.ui.control.Button
        preprocesssingButton          matlab.ui.control.Button
        readSACButton                 matlab.ui.control.Button
    end

    properties (Access = public)
        config           % 用于存放配置参数
        DataStruct       % 用于存放地震数据
        ccpStruct       % 用于存放ccp结果
        configEditStatus  % 'OK','Cancel', 或其他状态
        preprocessingStatus = 'Cancel' % 'OK','Cancel', 或其他状态
        rankReductionStatus = 'Cancel'% 'OK','Cancel', 或其他状态
        hkStackingStatus = 'Cancel'; % 'OK','Cancel', 或其他状态
        CCPStackingStatus % 'OK','Cancel', 或其他状态
    end

    properties (Access = private)
        seisout          % 叠加结果
        depth0           % 叠加结果对应的深度或时间坐标
        mohoStruct       % 存储莫霍面信息（如 stackCommonStationGather 返回时）
        Events           % 存放所有事件信息的结构数组
        Stations         % 存放所有台站信息的结构数组
        selectedEvid     % 存放选取的事件序号
        selectedStation  % 存放选取的台站序号
    end
    
    methods (Access = private)
        function tableData=getEventTable(app)
            % 1) 获取事件信息
            app.Events = getEvents(app.DataStruct); 
        
            % 2) 计算每个事件对应的“记录到的台站个数”
            %    可以先拿到 gather，再统计其长度
            nEvents = numel(app.Events);
            stationCounts = zeros(nEvents,1);
            for i = 1:nEvents
                evid   = app.Events(i).evid;
                gather = getCommonEventGather(app.DataStruct, evid);
                stationCounts(i) = numel(gather);
            end
        
            % 3) 构建要显示在 UITable 上的 Data（一般是 cell 或 table）
            % 假设想显示: [编号, lat, lon, depth, magnitude, stationCount]
            tableData = cell(nEvents,6);
            for i = 1:nEvents
                tableData{i,1} = app.Events(i).evid;                       % 事件ID (字符串)
                tableData{i,2} = sprintf('%.4f', app.Events(i).evla);      % 4位小数
                tableData{i,3} = sprintf('%.4f', app.Events(i).evlo);      % 4位小数
                tableData{i,4} = sprintf('%.2f', app.Events(i).evdp);      % 2位小数
                tableData{i,5} = sprintf('%.2f', app.Events(i).mag);       % 2位小数
                tableData{i,6} = sprintf('%d',   stationCounts(i));        % 整数
            end
        end

        function tableData=getStationTable(app)
            % 1) 获取事件信息
            app.Stations = getStations(app.DataStruct); 
        
            % 2) 计算每个台站对应的"记录到的事件个数"
            nStations = numel(app.Stations);
            eventCounts = zeros(nStations,1);
            for i = 1:nStations
                sta   = app.Stations(i).sta;
                gather = getCommonStationGather(app.DataStruct, sta);
                eventCounts(i) = numel(gather);
            end
        
            % 3) 构建要显示在 UITable 上的 Data（一般是 cell 或 table）
            % [台网, 台站，lat, lon, elevation, eventCount]
            tableData = cell(nStations,6);
            for i = 1:nStations
                tableData{i,1} = app.Stations(i).network;
                tableData{i,2} = app.Stations(i).sta;
                tableData{i,3} = sprintf('%.4f', app.Stations(i).stla);      % 4位小数
                tableData{i,4} = sprintf('%.4f', app.Stations(i).stlo);      % 4位小数
                tableData{i,5} = sprintf('%.2f', app.Stations(i).stel);      % 2位小数
                tableData{i,6} = sprintf('%d',   eventCounts(i));            % 整数
            end
        end


       function plotStationsInApp(app, selectedEvid, selectedStation)
           % PLOTSTATIONSINAPP  使用 GeographicAxes 同时绘制：
           %   - 所有台站 (灰色)
           %   - 选中事件台站 (红色)
           %   - 选中的单个台站 (蓝色)
           %
           % Inputs:
           %   app               - 主 App 对象
           %   selectedEvid      - 事件ID，用于高亮显示该事件的台站 (可为空)
           %   selectedStation   - 字符串，台站名，用于从台站表格点选后高亮 (可为空)
           %
           % 依赖:
           %   getStations(DataStruct)
           %   getCommonEventGather(DataStruct, evid)
           %   geoaxes, geoscatter, geolimits, geobasemap

           %------------------- 1. 获取 "所有台站" -------------------%
           allSt = getStations(app.DataStruct);  % struct array
           allLats  = [allSt.stla];
           allLons  = [allSt.stlo];
           allNames = {allSt.sta};

           %------------------- 2. 获取 "事件台站" -------------------%
           if ~exist('selectedEvid','var') || isempty(selectedEvid)
               eventLats  = [];
               eventLons  = [];
               eventNames = {};
           else
               gather  = getCommonEventGather(app.DataStruct, selectedEvid);
               eventSt = getStations(gather);
               eventLats  = [eventSt.stla];
               eventLons  = [eventSt.stlo];
               eventNames = {eventSt.sta};
           end


           %------------------- 3. 确定地理范围 (latlim, lonlim) -------------------%
           if ~isempty(allLats)
               latMin = min(allLats) - 0.25;
               latMax = max(allLats) + 0.25;
               lonMin = min(allLons) - 0.25;
               lonMax = max(allLons) + 0.25;
               if latMin == latMax, latMax = latMax + 0.5; end
               if lonMin == lonMax, lonMax = lonMax + 0.5; end
           else
               latMin = 30; latMax = 40;
               lonMin = 90; lonMax = 100;
           end

           %------------------- 4. 清空 Panel 并创建 geoaxes -------------------%
           delete(allchild(app.StationMapPanel));
           gax = geoaxes('Parent', app.StationMapPanel, ...
               'Units','normalized',...
               'Position',[0 0 1 1]);
           hold(gax, 'on');

           geobasemap(gax, 'topographic');
           geolimits(gax, [latMin latMax], [lonMin lonMax]);

           %------------------- 5. 绘制 "所有台站" (灰色) -------------------%
           if ~isempty(allLats)
               geoscatter(gax, allLats, allLons, ...
                   80, [0.7 0.7 0.7], '^','filled', ...  % 灰色
                   'MarkerEdgeColor','k');

%                % 标注台站名 (可酌情去掉以免拥挤)
%                for i = 1:numel(allLats)
%                    text(gax, allLats(i), allLons(i), allNames{i}, ...
%                        'FontSize',10, 'Color','k', ...
%                        'VerticalAlignment','bottom','HorizontalAlignment','left');
%                end
               
           end

           %------------------- 6. 高亮 "事件台站" (红色) -------------------%
           if ~isempty(eventLats)
               geoscatter(gax, eventLats, eventLons, ...
                   100, 'r','^','filled', ...
                   'MarkerEdgeColor','k');
%                for i = 1:numel(eventLats)
%                    text(gax, eventLats(i), eventLons(i), eventNames{i}, ...
%                        'FontSize',12, 'Color','r', ...
%                        'VerticalAlignment','bottom','HorizontalAlignment','right');
%                end
           end

           %------------------- 7. 若指定了 selectedStation, 再高亮选中台站 (蓝色) -------------------%
           if exist('selectedStation','var') && ~isempty(selectedStation)
               % 在 allSt 中查找站名
               idxSel = find(strcmp(allNames, selectedStation), 1);
               if ~isempty(idxSel)
                   selLat = allLats(idxSel);
                   selLon = allLons(idxSel);

                   % 以更醒目的符号/颜色绘制
                   geoscatter(gax, selLat, selLon, ...
                       120, 'b','^','filled', ...
                       'MarkerEdgeColor','k');

                   text(gax, selLat, selLon, selectedStation, ...
                       'FontSize',12, 'Color','b',...
                       'VerticalAlignment','bottom','HorizontalAlignment','center');
               else
                   % 如果没找到该 station, 也可给个提示
                   disp(['Station ', selectedStation, ' not found in allSt.']);
               end
           end

           %------------------- 8. 绘制断层 (可选) -------------------%
           faultDir = './visualization/faults/';
           if exist(faultDir,'dir')
               fault_files = dir([faultDir,'*txt']);
               for k = 1:length(fault_files)
                   faults = read_faults_gmt(fullfile(fault_files(k).folder, fault_files(k).name));
                   for l = 1:length(faults)
                       fault = faults{l};  % Nx2, [lon, lat]
                       geoplot(gax, fault(:,2), fault(:,1), 'LineWidth',2,'Color','k');
                   end
               end
           end
           hold(gax, 'off');
       end



        function plotAllEventsInApp(app,selectedEvid)
        % PLOTALLEVENTSINAPP 使用 GeographicAxes 在 app.EventMapPanel 中显示事件分布
        %
        % 假设:
        %   - app.Events 里有字段 evla, evlo, evid 等事件信息
        %   - app.EventMapPanel 是在 App Designer 中拖拽的一个 Panel，用于嵌入地图
            if nargin < 2
                selectedEvid = 'All';
            end
                
            % 1) 如果没有事件数据, 可以提前返回
            if isempty(app.Events)
                disp('No events to plot.');
                return;
            end
        
            % 2) 提取事件的纬经度与 ID
            lats = [app.Events.evla];
            lons = [app.Events.evlo];
            evidList = {app.Events.evid};
     
            % 3) 清空 MapPanel 的内容
            delete(allchild(app.EventMapPanel));
        
            % 4) 在 MapPanel 中创建一个 GeographicAxes
            gax = geoaxes('Parent', app.EventMapPanel, ...
                         'Units','normalized', ...
                         'Position',[0 0 1 1]);
            
            %（可选）如果您想在 GeographicAxes 上加交互，也可加:
            % GA.ButtonDownFcn = @(~,~) disp('Clicked on geoaxes');
        
            % 5) 绘制事件散点
            %    这里简单用 黑色圆点。MarkerSize 可视需求调整
            geoscatter(gax, lats, lons, 50, 'k', 'filled');
            hold(gax, 'on');

            % 绘制台站中心
%             Stations = getStations(app.DataStruct);
%             slat = mean([Stations.stla]);
%             slon = mean([Stations.stlo]);
            slat  = mean([app.Stations.stla]);
            slon  = mean([app.Stations.stlo]);
            geoscatter(gax, slat, slon, 100, 'r', '^','filled','MarkerEdgeColor','k');

            % 高亮选中的事件
            %    找到 selectedEvid 对应的索引
            idx = find(strcmp({app.Events.evid}, selectedEvid), 1);
            if ~isempty(idx)
                % 在图上用红色圈出
                geoscatter(gax, app.Events(idx).evla, app.Events(idx).evlo,50, ...
                    'ro','LineWidth',2);
               text(gax, app.Events(idx).evla, app.Events(idx).evlo, app.Events(idx).evid, ...
                   'FontSize',12, 'Color','r', ...
                   'VerticalAlignment','bottom','HorizontalAlignment','center');
            end
              
            % 6) 设置地图背景 (basemap)
            %    可选项有 'streets', 'satellite', 'topographic', 'darkwater' 等
            geobasemap(gax, 'topographic');
        
            % 7) 限制显示范围 (可选)
            %    如果您只想显示一定范围，可用 geolimits
            %    例如, geolimits(GA, [minLat maxLat], [minLon maxLon])
            % 这里示例：自动扩展到所有事件范围再加一点边界
            latMargin = 5;
            lonMargin = 5;
            latMin = min([slat lats]) - latMargin;
            latMax = max([slat lats]) + latMargin;
            lonMin = min([slon lons]) - lonMargin;
            lonMax = max([slon lons]) + lonMargin;
            geolimits(gax, [latMin latMax], [lonMin lonMax]);
        
            % 8) 设置标题
            title(gax, 'All Events (GeographicAxes)');
        
            % 9) 标注事件ID (可选)
%             for i = 1:numel(lats)
%                 text(gax, lats(i), lons(i), ['  ' evidList{i}], ...
%                     'FontSize',8, 'Color','k');
%             end
            
            hold(gax, 'off');
        end


        
        function updateEventStationTable(app)
            % 显示事件信息
            eventTableData = getEventTable(app);
            app.EventTable.Data = eventTableData;
            app.EventTable.ColumnName = {'Event ID','Latitude','Longitude','Depth','Mag','Station Count'};

            % 显示台站信息
            stationTableData = getStationTable(app);
            app.StationTable.Data = stationTableData;
            app.StationTable.ColumnName = {'Network','Station','Latitude','Longitude','Elevation','Event Count'};

        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: readSACButton
        function readSACButtonPushed(app, event)
%             % 0. 设置路径和参数
%             if isempty(app.config)
%                 uialert(app.UIFigure, '请先加载配置文件！', '错误','Icon','error');
%                 return;
%             end
%             try
%                 dataFolder = app.config.dataFolder;
%                 app.DataStruct = read_SAC(dataFolder);
%                 uialert(app.UIFigure, '数据读取完毕！', '提示');
%             catch ME
%                 uialert(app.UIFigure, ['数据读取失败: ', ME.message], '错误','Icon','error');
%             end
%             app.statusLabel.Text = '数据加载完成';
%             % 显示事件信息
%             loadEventList(app)
%             % 绘制事件位置
%             plotAllEventsInApp(app)
%             % 绘制台站位置
%             plotStationsInApp(app)
      
            if isempty(app.config)
                uialert(app.UIFigure, '请先加载配置文件！', '错误','Icon','error');
                return;
            end
            try
                % 1) 创建一个进度条对话框 (uiprogressdlg)
                d = uiprogressdlg(app.UIFigure, ...
                    'Title','读取数据', ...
                    'Message','正在读取 SAC 文件...', ...
                    'Cancelable','off', ...
                    'Value',0, ...
                    'Indeterminate','off');  
                
                % 2) read_SAC() 可能有多个文件循环，假设能传回或预先获取文件列表
                dataFolder = app.config.dataFolder;
                sacFiles = dir(fullfile(dataFolder, '*Z.*.SAC')); 
                nFiles = numel(sacFiles);
                DataStruct = [];
                for i = 1:nFiles
                    % 更新进度
                    d.Value = (i-1)/nFiles;
                    d.Message = sprintf('正在读取：%s (%d/%d)', sacFiles(i).name, i, nFiles);
                    drawnow;  % 强制界面刷新
        
                    % 调用 read_SAC_single(...) 或您自己的部分逻辑
                    sacFilePath = fullfile(dataFolder, sacFiles(i).name);
                    tmpStruct = read_SAC_single(sacFilePath);
                    if isempty(fieldnames(tmpStruct))
                        continue
                    end
                    DataStruct = [DataStruct tmpStruct];
                end
                
                % 最后
                d.Value = 1;
                d.Message = '所有 SAC 文件读取完成！';
                drawnow;
        
                pause(0.5);  % 给用户一点时间看进度到100%
                close(d);
        
                app.DataStruct = DataStruct;
                uialert(app.UIFigure, '数据读取完毕！', '提示');
            catch ME
                % 出错时记得关闭对话框
                if exist('d','var') && isvalid(d)
                    close(d);
                end
                uialert(app.UIFigure, ['数据读取失败: ', ME.message], '错误','Icon','error');
            end
            app.statusLabel.Text = '数据加载完成';

            % 更新事件台站列表
            updateEventStationTable(app)
%             % 显示事件信息
%             eventTableData = getEventTable(app);
%             app.EventTable.Data = eventTableData;
%             app.EventTable.ColumnName = {'Event ID','Latitude','Longitude','Depth','Mag','Station Count'};
% 
%             % 显示台站信息
%             stationTableData = getStationTable(app);
%             app.StationTable.Data = stationTableData;
%             app.StationTable.ColumnName = {'Network','Station','Latitude','Longitude','Elevation','Event Count'};
            % 绘制地震位置
            plotAllEventsInApp(app)
            % 绘制台站位置
            plotStationsInApp(app)
        end

        % Button pushed function: preprocesssingButton
        function preprocesssingButtonPushed(app, event)
%             % 2. 预处理
%             if isempty(app.DataStruct)
%                 uialert(app.UIFigure, '尚未读取数据！', '错误','Icon','error');
%                 return;
%             end
%             if ~isfield(app.config,'PreprocessingParam')
%                 uialert(app.UIFigure, '配置中缺少 PreprocessingParam！', '错误','Icon','error');
%                 return;
%             end
%             try
%                 app.DataStruct = preprocessing(app.DataStruct, app.config.PreprocessingParam);
%                 uialert(app.UIFigure, '预处理完成！', '提示');
%             catch ME
%                 uialert(app.UIFigure, ['预处理失败: ', ME.message], '错误','Icon','error');
%             end
%             app.statusLabel.Text = '预处理完成';


            % 当用户点击 Preprocessing 按键时执行
            try
                if isempty(app.DataStruct)
                    uialert(app.UIFigure,'尚未加载数据！','提示');
                    return;
                end

                % 如果 param 还在 app.config 里
                if ~isfield(app.config,'PreprocessingParam')
                    uialert(app.UIFigure,'未找到 PreprocessingParam!','提示');
                    return;
                end
                defaultParam = app.config.PreprocessingParam;

                % 创建并打开子App
                preprocessingApp = PreprocessingApp();
                preprocessingApp.mainApp = app;
                preprocessingApp.DataStruct = app.DataStruct;
                preprocessingApp.loadConfigData(defaultParam);
                % 让子App知道 mainApp 的引用 -> preprocessingApp.mainApp = app;

                % 如果要阻塞主App，可:
                waitfor(preprocessingApp.UIFigure);

                % 现在子App 关闭后，若它更新了 app.DataStruct，主App 即可用
                
                switch app.preprocessingStatus
                    case 'OK'
                        uialert(app.UIFigure,'Preprocessing 完成！','提示');
                        app.statusLabel.Text = 'Preprocessing 完成';
                        % 更新事件台站列表
                        updateEventStationTable(app)
                    case 'Cancel'
                        uialert(app.UIFigure,'Preprocessing 取消!','提示');
                        app.statusLabel.Text = 'Preprocessing 取消';
                    otherwise
                        uialert(app.UIFigure,'Known operation...','提示');
                end
                
            catch ME
                uialert(app.UIFigure, ...
                    ['打开Preprocessing参数设置时出现错误: ',ME.message], ...
                    '错误','Icon','error');
            end


%             if isempty(app.DataStruct)
%                 uialert(app.UIFigure, '尚未读取数据！', '错误','Icon','error');
%                 return;
%             end
%             if ~isfield(app.config,'PreprocessingParam')
%                 uialert(app.UIFigure, '配置中缺少 PreprocessingParam！', '错误','Icon','error');
%                 return;
%             end
%             try
%                 % 1) 创建一个进度条对话框 (uiprogressdlg)
%                 d = uiprogressdlg(app.UIFigure, ...
%                     'Title','预处理', ...
%                     'Message','正在进行预处理，请稍候...', ...
%                     'Cancelable','off', ...
%                     'Value',0, ...
%                     'Indeterminate','off');
% 
%                 nFiles = length(app.DataStruct);
%                 
%                 for i = 1:nFiles
%                     % 更新进度
%                     d.Value = (i-1)/nFiles;
%                     d.Message = sprintf('正在预处理：第 (%d/%d) 条数据', i, nFiles);
%                     drawnow;  % 强制界面刷新
%        
%                     % 调用 preprocessing_single(...)
%                     [app.DataStruct(i),removeFlag] = preprocessing_single(app.DataStruct(i), app.config.PreprocessingParam);
%                 end
%                 
%                 % 最后
%                 d.Value = 1;
%                 d.Message = '所有数据预处理完成！';
%                 drawnow;
%         
%                 pause(0.5);  % 给用户一点时间看进度到100%
%                 close(d);
%         
%                 uialert(app.UIFigure, '预处理完成！', '提示');
%             catch ME
%                 if exist('d','var') && isvalid(d)
%                     close(d);
%                 end
%                 uialert(app.UIFigure, ['预处理失败: ', ME.message], '错误','Icon','error');
%             end
%             app.statusLabel.Text = '预处理完成';

        end

        % Button pushed function: deconvButton
        function deconvButtonPushed(app, event)
            % 3. 计算接收函数
%             if isempty(app.DataStruct)
%                 uialert(app.UIFigure, '尚未读取并预处理数据！', '错误','Icon','error');
%                 return;
%             end
%             if ~isfield(app.config,'DeconvParam')
%                 uialert(app.UIFigure, '配置中缺少 DeconvParam！', '错误','Icon','error');
%                 return;
%             end
%             try
%                 app.DataStruct = deconv(app.DataStruct, app.config.DeconvParam);
%                 uialert(app.UIFigure, '接收函数计算完成！', '提示');
%             catch ME
%                 uialert(app.UIFigure, ['接收函数计算失败: ', ME.message], '错误','Icon','error');
%             end
% 
%             app.statusLabel.Text = '接收函数计算完成';

            if isempty(app.DataStruct)
                uialert(app.UIFigure, '尚未读取并预处理数据！', '错误','Icon','error');
                return;
            end
            if ~isfield(app.config,'DeconvParam')
                uialert(app.UIFigure, '配置中缺少 DeconvParam！', '错误','Icon','error');
                return;
            end
            try
                d = uiprogressdlg(app.UIFigure,...
                    'Title','计算接收函数',...
                    'Message','正在进行去卷积处理...',...
                    'Value',0, ...
                    'Indeterminate','off',...
                    'Cancelable','off');
                nFiles = length(app.DataStruct);
                noresults = [];
                for i = 1:nFiles
                    % 更新进度
                    d.Value = (i-1)/nFiles;
                    d.Message = sprintf('正在反褶积：第 (%d/%d) 条数据', i, nFiles);
                    drawnow;  % 强制界面刷新
       
                    % 调用 deconv_single(...)
                    [app.DataStruct(i),noresult] = deconv_single(app.DataStruct(i), app.config.DeconvParam);
                    if noresult
                        noresults = [noresults,i];
                    end
                end
                
                % 最后
                d.Value = 1;
                d.Message = '接收函数计算完成！';
                drawnow;


                pause(0.5);  % 给用户一点时间看进度到100%
                close(d);

                uialert(app.UIFigure, '接收函数计算完成！', '提示');
            catch ME
                if exist('d','var') && isvalid(d)
                    close(d);
                end
                uialert(app.UIFigure, ['接收函数计算失败: ', ME.message], '错误','Icon','error');
            end

            if ~isempty(noresults)
                % If desired, remove the invalid or unsuccessful records
                app.DataStruct(noresults) = [];
                fprintf('Removed %d traces with no valid decon results.\n', length(noresult));
            end

            app.statusLabel.Text = '接收函数计算完成';
            % 更新事件台站列表
            updateEventStationTable(app)
%             % 更新显示事件信息列表
%             eventTableData = getEventTable(app);
%             app.EventTable.Data = eventTableData;
%             app.EventTable.ColumnName = {'Event ID','Latitude','Longitude','Depth','Mag','Station Count'};
%             % 更新显示台站信息列表
%             stationTableData = getStationTable(app);
%             app.StationTable.Data = stationTableData;
%             app.StationTable.ColumnName = {'Network','Station','Latitude','Longitude','Elevation','Event Count'};
        end

        % Button pushed function: stackButton
        function stackButtonPushed(app, event)
            % 4. 叠加接收函数
            if isempty(app.DataStruct)
                uialert(app.UIFigure, '尚未进行接收函数计算！', '错误','Icon','error');
                return;
            end
            try
                [app.seisout, app.depth0, app.mohoStruct] = stackCommonStationGather(app.DataStruct);
                uialert(app.UIFigure, '接收函数叠加完成！', '提示');
            catch ME
                uialert(app.UIFigure, ['叠加失败: ', ME.message], '错误','Icon','error');
            end
            app.statusLabel.Text = '接收函数叠加完成';
        end

        % Button pushed function: setupPathsButton
        function setupPathsButtonPushed(app, event)
            try
                setupPaths();
                uialert(app.UIFigure, '项目路径已设置完毕！', '提示');
            catch ME
                uialert(app.UIFigure, ['设置路径失败: ', ME.message], '错误', 'Icon','error');
            end
            app.statusLabel.Text = '路径设置完成';
        end

        % Button pushed function: loadConfigButton
        function loadConfigButtonPushed(app, event)
            try
                defaultCfg = loadConfig();
                app.configEditStatus = 'Initial';  % 初始值
            
                configEditor = ConfigEditorApp();
                configEditor.mainApp = app;
                configEditor.loadConfigData(defaultCfg);
            
                waitfor(configEditor.UIFigure);
            
                switch app.configEditStatus
                    case 'OK'
                        uialert(app.UIFigure,'配置已加载并更新！','提示');
                        app.statusLabel.Text = '配置已加载并更新';
                    case 'Cancel'
                        uialert(app.UIFigure,'已取消更新配置，保留默认设置','提示');
                        app.statusLabel.Text = '已取消更新配置，保留默认设置';
                    otherwise
                        uialert(app.UIFigure,'未知操作...','提示');
                end
            
            catch ME
                uialert(app.UIFigure, ...
                    ['加载或编辑配置时出现错误: ', ME.message], ...
                    '错误','Icon','error');
            end

        end

        % Button pushed function: plotWaveformsButton
        function plotWaveformsButtonPushed(app, event)
            if isempty(app.DataStruct)
                uialert(app.UIFigure, '尚未有可用数据！', '错误','Icon','error');
                return;
            end
            answer = inputdlg({'请输入要绘制的道索引：'}, ...
                              '绘制单条记录', ...
                              [1 35], ...
                              {'100'});
            if isempty(answer)
                return; % 用户取消
            end
            traceIndex = str2double(answer{1});
            if isnan(traceIndex) || traceIndex <= 0
                uialert(app.UIFigure, '道索引必须为正整数！', '错误','Icon','error');
                return;
            end
            try
%                 figure('Name','Single Trace','NumberTitle','off');
                plotWaveforms(app.DataStruct, traceIndex);
            catch ME
                uialert(app.UIFigure, ['绘图失败: ', ME.message], '错误','Icon','error');
            end
        end

        % Button pushed function: plotComonStationGatherButton
        function plotComonStationGatherButtonPushed(app, event)
            if isempty(app.DataStruct)
                uialert(app.UIFigure, '尚未有可用数据！', '错误','Icon','error');
                return;
            end
            try
                stations = getStations(app.DataStruct);
                stationList = {stations.sta};
                [selIndex, ok] = listdlg('PromptString','请选择台站：',...
                                         'SelectionMode','single',...
                                         'ListString',stationList);
                if ~ok, return; end
                chosenStation = stationList{selIndex};
                
%                 figure('Name','Common Station Gather','NumberTitle','off');
                plotCommonStationGather(app.DataStruct, chosenStation);
            catch ME
                uialert(app.UIFigure, ['绘图失败: ', ME.message], '错误','Icon','error');
            end
        end

        % Button pushed function: plotCommonEventGatherButton
        function plotCommonEventGatherButtonPushed(app, event)
        % 该函数演示根据用户选择的事件编号（或 "All"）绘制事件道集。
        % 支持一次性绘制多个事件或所有事件。
    
        if isempty(app.DataStruct)
            uialert(app.UIFigure, '尚未有可用数据！', '错误','Icon','error');
            return;
        end
        
        % 1) 获取事件列表
        try
            events = getEvents(app.DataStruct);
            eventList = {events.evid};  % cell array of strings
        catch ME
            uialert(app.UIFigure, ['无法获取事件列表: ', ME.message], '错误','Icon','error');
            return;
        end

        % 2) 构造带有 "All" 的选项列表
        extendedEventList = [{'All'}, eventList];

        % 3) 弹出多选列表对话框
        [selIndex, ok] = listdlg('PromptString',   '请选择事件：',...
            'ListString',     extendedEventList,...
            'SelectionMode',  'multiple',...
            'ListSize',       [200 300]);  % 可自行调整大小
        if ~ok
            return;  % 用户点击“取消”
        end

        % 4) 选择一个文件夹保存绘图结果（如果需要）
%         outputDir = uigetdir(pwd, '选择输出图像保存文件夹');
%         if outputDir == 0
            outputDir = './figures'; % 用户取消选择目录
%         end

        % 5) 判断用户是否选中了 "All"
        if any(selIndex == 1)
            % 用户选中了 "All"
            % 等效于在列表中选中所有实际的事件 ID
            selectedEventIDs = eventList;  % 全部
        else
            % 用户只选中了具体的事件（单个或多个）
            % 索引需要减去1，因为第1项是 "All"
            realIndex = selIndex - 1;
            selectedEventIDs = eventList(realIndex);
        end

        % 6) 绘制并保存
        try
            for iE = 1:numel(selectedEventIDs)
                evid = selectedEventIDs{iE};

                % 创建一个不可见的 figure 以便批量输出
%                 f = figure('Visible','off','NumberTitle','off',...
%                     'Name',['Event ', evid]);

                % 根据需要指定绘图模式，如 'trace'、'wiggle' 等
                plotCommonEventGather(app.DataStruct, evid, 'trace','wigb');

                % 保存图像文件 (若安装 export_fig 也可使用 export_fig)
                filename = fullfile(outputDir, [evid, '.png']);
                saveas(gcf, filename);

                close(gcf);
            end

            uialert(app.UIFigure, '事件道集绘制并保存完成！','提示');
        catch ME
            uialert(app.UIFigure, ['绘图或保存失败: ', ME.message], '错误','Icon','error');
        end
        end

        % Button pushed function: plotStationsEventsButton
        function plotStationsEventsButtonPushed(app, event)
            if isempty(app.DataStruct)
                uialert(app.UIFigure, '尚未有可用数据！', '错误','Icon','error');
                return;
            end
            try
                demFile = 'Qaidam_DEM.mat';
                if ~exist(demFile,'file')
                    warndlg(['未找到 DEM 文件: ', demFile, ...
                        '。仅绘制台站和事件位置。'],'警告');
                    demFile = [];
                end

                % 绘制台站分布
                plotStations(app.DataStruct, demFile);

                pause(0.05);
                
                % 绘制事件分布
                plotEvents(app.DataStruct);

            catch ME
                uialert(app.UIFigure, ['绘图失败: ', ME.message], '错误','Icon','error');
            end
        end

        % Cell selection callback: EventTable
        function EventTableCellSelection(app, event)
            indices = event.Indices;
            % Value changed function: EventsTable
            if isempty(indices)
                return; % 用户可能取消选中
            end
            selectedRow = indices(1);

            % 获取该行对应的事件信息
            app.selectedEvid = app.Events(selectedRow).evid;
            % 进行绘图
            plotAllEventsInApp(app,app.selectedEvid)
            plotStationsInApp(app,app.selectedEvid);
        end

        % Selection change function: InfoTabGroup
        function InfoTabGroupSelectionChanged(app, event)
            selectedTab = app.InfoTabGroup.SelectedTab;
            switch selectedTab.Title
                case 'Event Table'
                    % 显示地震事件表格
                    app.EventTable.Visible = 'on';
                    app.StationTable.Visible = 'off';

                case 'Station Table'
                    % 显示台站表格
                    app.StationTable.Visible = 'on';
                    app.EventTable.Visible = 'off';  
            end
        end

        % Cell selection callback: StationTable
        function StationTableCellSelection(app, event)
            idx = event.Indices;
            if isempty(idx)
                return;
            end
            rowIdx = idx(1);  % 选中行
            
            % 假设台站名称在第2列(Station)
            app.selectedStation = app.StationTable.Data{rowIdx,2};
        
            % 在这里调用 plotStationsInApp，传入 selectedStation
            % 也可以决定要不要传 selectedEvid
            plotStationsInApp(app, app.selectedEvid, app.selectedStation);
        end

        % Menu selected function: LoadMenu
        function onLoadMenuSelected(app, event)
            try
                % 弹出“打开文件”对话框，让用户选择数据文件
                [filename, pathname] = uigetfile({'*.mat;*.csv;*.txt','Data Files'}, ...
                    '请选择要加载的数据文件');
                if isequal(filename,0) || isequal(pathname,0)
                    % 用户取消选择
                    return;
                end

                fullpath = fullfile(pathname, filename);
                % 根据文件后缀或您项目需求来读入数据
                % 这里只示范 .mat 的情况
                [~,~,ext] = fileparts(filename);
                switch lower(ext)
                    case '.mat'
                        loaded = load(fullpath);  % 载入 .mat
                        if isfield(loaded,'DataStruct')
                            app.DataStruct = loaded.DataStruct;
                        else
                            % 如果 mat 里不是这个字段，可自行处理
                            app.DataStruct = loaded;
                        end
                    otherwise
                        % 如 .csv 等, 可以用 readtable, csvread, textread, etc.
                        % 简单示例:
                        if strcmp(ext,'.csv')
                            T = readtable(fullpath);
                            app.DataStruct = table2struct(T);
                        elseif strcmp(ext,'.txt')
                            % ...
                        end
                end

                % 根据加载结果做相应处理（如刷新界面、更新表格等）
                uialert(app.UIFigure, '数据加载成功！', '提示');

                % 显示事件信息
                eventTableData = getEventTable(app);
                app.EventTable.Data = eventTableData;
                app.EventTable.ColumnName = {'Event ID','Latitude','Longitude','Depth','Mag','Station Count'};
    
                % 显示台站信息
                stationTableData = getStationTable(app);
                app.StationTable.Data = stationTableData;
                app.StationTable.ColumnName = {'Network','Station','Latitude','Longitude','Elevation','Event Count'};
                % 绘制地震位置
                plotAllEventsInApp(app)
                % 绘制台站位置
                plotStationsInApp(app)
            catch ME
                % 出错处理
                uialert(app.UIFigure, ...
                    ['加载数据时出现错误: ', ME.message], ...
                    '错误','Icon','error');
            end
        end

        % Menu selected function: SaveMenu
        function onSaveMenuSelected(app, event)
            if isempty(app.DataStruct)
                uialert(app.UIFigure, '暂无数据可保存！', '提示');
                return;
            end

            try
                [filename, pathname] = uiputfile('*.mat','选择保存文件');
                if isequal(filename,0) || isequal(pathname,0)
                    % 用户取消
                    return;
                end

                fullpath = fullfile(pathname, filename);

                % 简单示例：把 app.DataStruct 存为 .mat
                DataStruct=app.DataStruct;
                save(fullpath,'DataStruct','-v7.3');

                uialert(app.UIFigure, ['数据已保存到: ', fullpath], '提示');

            catch ME
                uialert(app.UIFigure, ...
                    ['保存数据时出现错误: ', ME.message], ...
                    '错误','Icon','error');
            end
        end

        % Menu selected function: RankReductionMenu
        function onRankReductionMenuSelected(app, event)
            % 当用户点击 Tools->Rank Reduction 菜单时执行
            try
                % 如果您希望在这里就直接对 app.DataStruct 做 rankReduction：
                %   1) 需要先拿到 param
                %   2) rankReduction(app.DataStruct, param)
                %   3) 更新 app.DataStruct

                % 但是通常您需要让用户先编辑参数 -> 故建议弹出一个 "RankReductionApp"
                % 传入 app.DataStruct, app.config.RankReductionParam, etc.

                if isempty(app.DataStruct)
                    uialert(app.UIFigure,'尚未加载数据！','提示');
                    return;
                end

                % 如果 param 还在 app.config 里
                if ~isfield(app.config,'RankReductionParam')
                    uialert(app.UIFigure,'未找到 RankReductionParam!','提示');
                    return;
                end
                defaultParam = app.config.RankReductionParam;

                % 创建并打开子App
                rankApp = RankReductionApp();
                rankApp.mainApp = app;
                rankApp.DataStruct = app.DataStruct;
                rankApp.loadConfigData(defaultParam);
                % 让子App知道 mainApp 的引用 -> rankApp.mainApp = app;

                % 如果要阻塞主App，可:
                waitfor(rankApp.UIFigure);

                % 现在子App 关闭后，若它更新了 app.DataStruct，主App 即可用
                % 也可做提示
%                 uialert(app.UIFigure,'Rank Reduction 处理完成或已取消！','提示');
                
                switch app.rankReductionStatus
                    case 'OK'
                        uialert(app.UIFigure,'Rank reduction 处理完成！','提示');
                        app.statusLabel.Text = 'Rank reduction 处理完成';
                    case 'Cancel'
                        uialert(app.UIFigure,'Rank reduction 取消!','提示');
                        app.statusLabel.Text = 'Rank reduction 取消';
                    otherwise
                        uialert(app.UIFigure,'Known operation...','提示');
                end
                % 更新事件台站列表
                updateEventStationTable(app)

            catch ME
                uialert(app.UIFigure, ...
                    ['打开Rank Reduction参数设置时出现错误: ',ME.message], ...
                    '错误','Icon','error');
            end

           
        end


        % Menu selected function: HkStackingMenu
        function HkStackingMenuSelected(app, event)
            try
                % 检查是否加载了数据
                if isempty(app.DataStruct)
                    uialert(app.UIFigure,'尚未加载数据！','提示');
                    return;
                end
                
                % 如果你有 HK stacking 的参数配置，并存放在 app.config 中，可取出默认参数
                if isfield(app.config, 'HKStackingParam')
                    defaultParam = app.config.HKStackingParam;
                else
                    defaultParam = [];
                end
        
                % 创建并打开子 App
                hkApp = HkStackingApp();     % 调用你实现的 HKStackingApp
                hkApp.mainApp = app;          % 传入主 App 的引用
                hkApp.DataStruct = app.DataStruct;  % 传入数据
%                 hkApp.refreshUI(); 
                % 如果 HKStackingApp 实现了加载参数的函数，则调用它
                if ~isempty(defaultParam) && ismethod(hkApp, 'loadConfigData')
                    hkApp.loadConfigData(defaultParam);
                end
                
                % 如果需要让主 App 阻塞，等待子 App 关闭，则调用 waitfor
                waitfor(hkApp.UIFigure);

                % 子 App 关闭后，根据 hkStackingStatus 更新主 App 状态
                switch app.hkStackingStatus
                        case 'OK'
                            uialert(app.UIFigure,'HK Stacking 处理完成！','提示');
                            app.statusLabel.Text = 'HK Stacking 处理完成';
                        case 'Cancel'
                            uialert(app.UIFigure,'HK Stacking 取消！','提示');
                            app.statusLabel.Text = 'HK Stacking 取消';
                        otherwise
                            uialert(app.UIFigure,'未知操作...','提示');
                end
            catch ME
                uialert(app.UIFigure, ...
                    ['打开HK Stacking参数设置时出现错误: ', ME.message], ...
                    '错误','Icon','error');
            end

        % Menu selected function: CCPMenu
        function CCPMenuSelected(app, event)
             % 当用户点击 imaging->CCP 菜单时执行
            try

                % 创建并打开子App
                ccpApp = CCPStackingApp();
                ccpApp.mainApp = app;
                ccpApp.DataStruct = app.DataStruct;
                % ccpApp.loadConfigData(defaultParam);
                % 让子App知道 mainApp 的引用 -> rankApp.mainApp = app;

                % 如果要阻塞主App，可:
                waitfor(ccpApp.UIFigure);

                switch app.CCPStackingStatus
                    case 'OK'
                        uialert(app.UIFigure,'CCP stacking 处理完成！','提示');
                        app.statusLabel.Text = 'CCP stacking 处理完成';
                    case 'Cancel'
                        uialert(app.UIFigure,'CCP stacking 取消!','提示');
                        app.statusLabel.Text = 'CCP stacking 取消';
                    otherwise
                        uialert(app.UIFigure,'Known operation...','提示');
                end        

            catch ME
                uialert(app.UIFigure, ...
                    ['打开CommonConversionPoint参数设置时出现错误: ',ME.message], ...
                    '错误','Icon','error');
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [198 97 1040 650];
            app.UIFigure.Name = 'MATLAB App';

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadMenu
            app.LoadMenu = uimenu(app.FileMenu);
            app.LoadMenu.MenuSelectedFcn = createCallbackFcn(app, @onLoadMenuSelected, true);
            app.LoadMenu.Text = '加载 (Load)';

            % Create SaveMenu
            app.SaveMenu = uimenu(app.FileMenu);
            app.SaveMenu.MenuSelectedFcn = createCallbackFcn(app, @onSaveMenuSelected, true);
            app.SaveMenu.Text = '保存 (Save)';

            % Create ToolsMenu
            app.ToolsMenu = uimenu(app.UIFigure);
            app.ToolsMenu.Text = 'Tools';

            % Create RadonTransformMenu
            app.RadonTransformMenu = uimenu(app.ToolsMenu);
            app.RadonTransformMenu.Text = 'Radon Transform';

            % Create RankReductionMenu
            app.RankReductionMenu = uimenu(app.ToolsMenu);
            app.RankReductionMenu.MenuSelectedFcn = createCallbackFcn(app, @onRankReductionMenuSelected, true);
            app.RankReductionMenu.Text = 'Rank Reduction';

            % Create ImagingMenu
            app.ImagingMenu = uimenu(app.UIFigure);
            app.ImagingMenu.Text = 'Imaging';

            % Create HkStackingMenu
            app.HkStackingMenu = uimenu(app.ImagingMenu);
            app.HkStackingMenu.MenuSelectedFcn = createCallbackFcn(app, @HkStackingMenuSelected, true);
            app.HkStackingMenu.Text = 'H-k Stacking';

            % Create CCPMenu
            app.CCPMenu = uimenu(app.ImagingMenu);
            app.CCPMenu.MenuSelectedFcn = createCallbackFcn(app, @CCPMenuSelected, true);
            app.CCPMenu.Text = 'CCP';

            % Create Migration2DMenu
            app.Migration2DMenu = uimenu(app.ImagingMenu);
            app.Migration2DMenu.Text = 'Migration 2D';

            % Create Migration3DMenu
            app.Migration3DMenu = uimenu(app.ImagingMenu);
            app.Migration3DMenu.Text = 'Migration 3D';

            % Create EditingMenu
            app.EditingMenu = uimenu(app.UIFigure);
            app.EditingMenu.Text = 'Editing';

            % Create ViewMenu
            app.ViewMenu = uimenu(app.UIFigure);
            app.ViewMenu.Text = 'View';

            % Create readSACButton
            app.readSACButton = uibutton(app.UIFigure, 'push');
            app.readSACButton.ButtonPushedFcn = createCallbackFcn(app, @readSACButtonPushed, true);
            app.readSACButton.Position = [17 493 160 25];
            app.readSACButton.Text = '读入数据（readSAC)';

            % Create preprocesssingButton
            app.preprocesssingButton = uibutton(app.UIFigure, 'push');
            app.preprocesssingButton.ButtonPushedFcn = createCallbackFcn(app, @preprocesssingButtonPushed, true);
            app.preprocesssingButton.Position = [17 455 160 25];
            app.preprocesssingButton.Text = '预处理（preprocesssing)';

            % Create deconvButton
            app.deconvButton = uibutton(app.UIFigure, 'push');
            app.deconvButton.ButtonPushedFcn = createCallbackFcn(app, @deconvButtonPushed, true);
            app.deconvButton.Position = [17 417 160 25];
            app.deconvButton.Text = '计算接收函数 (deconv)';

            % Create stackButton
            app.stackButton = uibutton(app.UIFigure, 'push');
            app.stackButton.ButtonPushedFcn = createCallbackFcn(app, @stackButtonPushed, true);
            app.stackButton.Position = [17 379 160 25];
            app.stackButton.Text = '叠加接收函数 (stack)';

            % Create plotWaveformsButton
            app.plotWaveformsButton = uibutton(app.UIFigure, 'push');
            app.plotWaveformsButton.ButtonPushedFcn = createCallbackFcn(app, @plotWaveformsButtonPushed, true);
            app.plotWaveformsButton.Position = [17 253 160 25];
            app.plotWaveformsButton.Text = '绘制单条记录';

            % Create visualizationOptionLabel
            app.visualizationOptionLabel = uilabel(app.UIFigure);
            app.visualizationOptionLabel.FontWeight = 'bold';
            app.visualizationOptionLabel.Position = [17 282 114 22];
            app.visualizationOptionLabel.Text = '可视化操作：';

            % Create plotComonStationGatherButton
            app.plotComonStationGatherButton = uibutton(app.UIFigure, 'push');
            app.plotComonStationGatherButton.ButtonPushedFcn = createCallbackFcn(app, @plotComonStationGatherButtonPushed, true);
            app.plotComonStationGatherButton.Position = [17 216 160 25];
            app.plotComonStationGatherButton.Text = '绘制台站道集';

            % Create plotCommonEventGatherButton
            app.plotCommonEventGatherButton = uibutton(app.UIFigure, 'push');
            app.plotCommonEventGatherButton.ButtonPushedFcn = createCallbackFcn(app, @plotCommonEventGatherButtonPushed, true);
            app.plotCommonEventGatherButton.Position = [17 179 160 25];
            app.plotCommonEventGatherButton.Text = '绘制事件道集';

            % Create plotStationsEventsButton
            app.plotStationsEventsButton = uibutton(app.UIFigure, 'push');
            app.plotStationsEventsButton.ButtonPushedFcn = createCallbackFcn(app, @plotStationsEventsButtonPushed, true);
            app.plotStationsEventsButton.Position = [17 142 160 25];
            app.plotStationsEventsButton.Text = '绘制台站&事件分布';

            % Create loadConfigButton
            app.loadConfigButton = uibutton(app.UIFigure, 'push');
            app.loadConfigButton.ButtonPushedFcn = createCallbackFcn(app, @loadConfigButtonPushed, true);
            app.loadConfigButton.Position = [17 532 160 25];
            app.loadConfigButton.Text = '加载配置（loadConfig)';

            % Create setupPathsButton
            app.setupPathsButton = uibutton(app.UIFigure, 'push');
            app.setupPathsButton.ButtonPushedFcn = createCallbackFcn(app, @setupPathsButtonPushed, true);
            app.setupPathsButton.Position = [17 571 160 25];
            app.setupPathsButton.Text = '设置路径（setupPaths)';

            % Create processingWorkflowLabel
            app.processingWorkflowLabel = uilabel(app.UIFigure);
            app.processingWorkflowLabel.FontWeight = 'bold';
            app.processingWorkflowLabel.Position = [17 607 114 22];
            app.processingWorkflowLabel.Text = '处理流程：';

            % Create statusLabel
            app.statusLabel = uilabel(app.UIFigure);
            app.statusLabel.BackgroundColor = [1 1 1];
            app.statusLabel.FontWeight = 'bold';
            app.statusLabel.FontColor = [0.149 0.149 0.149];
            app.statusLabel.Position = [84 607 205 22];
            app.statusLabel.Text = '';

            % Create EventMapPanel
            app.EventMapPanel = uipanel(app.UIFigure);
            app.EventMapPanel.Title = 'Event Map Panel';
            app.EventMapPanel.Position = [228 349 354 254];

            % Create StationMapPanel
            app.StationMapPanel = uipanel(app.UIFigure);
            app.StationMapPanel.Title = 'Station Map Panel';
            app.StationMapPanel.Position = [603 349 354 254];

            % Create ButtonGroup
            app.ButtonGroup = uibuttongroup(app.UIFigure);
            app.ButtonGroup.Title = 'Button Group';
            app.ButtonGroup.Position = [-29 701 2 2];

            % Create InfoTabGroup
            app.InfoTabGroup = uitabgroup(app.UIFigure);
            app.InfoTabGroup.SelectionChangedFcn = createCallbackFcn(app, @InfoTabGroupSelectionChanged, true);
            app.InfoTabGroup.Position = [228 56 728 283];

            % Create EventTab
            app.EventTab = uitab(app.InfoTabGroup);
            app.EventTab.Title = 'Event Table';

            % Create EventTable
            app.EventTable = uitable(app.EventTab);
            app.EventTable.ColumnName = {'Event ID'; 'Longitude'; 'Latitude'; 'Depth'; 'Magnitude'; 'Station Count'};
            app.EventTable.RowName = {};
            app.EventTable.ColumnEditable = true;
            app.EventTable.CellSelectionCallback = createCallbackFcn(app, @EventTableCellSelection, true);
            app.EventTable.Position = [1 1 727 252];

            % Create StationTab
            app.StationTab = uitab(app.InfoTabGroup);
            app.StationTab.Title = 'Station Table';

            % Create StationTable
            app.StationTable = uitable(app.StationTab);
            app.StationTable.ColumnName = {'Network'; 'Station'; 'Longitude'; 'Latitude'; 'Elevation'; 'Event Count'};
            app.StationTable.RowName = {};
            app.StationTable.ColumnEditable = true;
            app.StationTable.CellSelectionCallback = createCallbackFcn(app, @StationTableCellSelection, true);
            app.StationTable.Position = [1 1 727 252];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DenseArrayToolkit_GUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
