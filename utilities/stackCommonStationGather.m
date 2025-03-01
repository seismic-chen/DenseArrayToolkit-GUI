function [seisout, depth0, mohoStruct] = stackCommonStationGather(DataStruct)
    % stackCommonStationGather - 叠加台站接收函数，提取莫霍面深度并绘制结果
    %
    % 输入:
    %   DataStruct - 包含台站和事件数据的结构体
    %
    % 输出:
    %   seisout    - 叠加后的接收函数剖面 [nz x nsta]
    %   depth0     - 深度轴 [nz x 1]
    %   mohoStruct - 包含莫霍面深度信息的结构体

    % 1) 获取台站列表
    station = getStations(DataStruct); 
    stationList = {station.sta};
    nsta = length(stationList);

    % 2) 预分配数组
    [seisCell, pArr, stlas, stlos] = initializeArrays(nsta);

    % 3) 加载速度模型
    [z, ~, vp, vs] = loadVelocityModel();
    zmax = 100; dz = 0.5;  % 深度范围和采样间隔

    % 4) 遍历台站，叠加波形
    for n = 1:nsta
        gather = getCommonStationGather(DataStruct, stationList{n});
        if isempty(gather)
            continue;  % 跳过空数据
        end

        % 提取台站信息
        stlas(n) = gather(1).StationInfo.stla;
        stlos(n) = gather(1).StationInfo.stlo;

        % 叠加波形
        [seisCell{n}, pArr(n)] = stackWaveforms(gather, dz, zmax, z, vp, vs);
    end

    % 5) 转换为矩阵并平滑
    ngrid_x = 8;
    ngrid_y = 4;
%     seisout = smoothSeismicData(seisCell,ngrid_x,ngrid_y);
    seisout = cell2mat(seisCell);
    % 6) 提取深度轴
    depth0 = (0:dz:zmax)';  % 深度轴 [nz x 1]

    % 7) 提取莫霍面深度
    [dmoho, amoho] = extractMohoDepth(seisout, depth0);

    % 8) 剔除离群值并插值
    [Fmoho, idx] = interpolateMohoDepth(stlas, stlos, dmoho);

    % 9) 绘制莫霍面深度图
    plotMohoMap(stlas, stlos, dmoho, Fmoho);

    % 10) 绘制叠加剖面
    plotStackedSection(seisout, depth0, dmoho);

    % 11) 输出结果
    mohoStruct = createMohoStruct(stlas, stlos, dmoho, amoho);
end

%% 子函数：初始化数组
function [seisCell, pArr, stlas, stlos] = initializeArrays(nsta)
    seisCell = cell(1, nsta);  % 存储叠加后的接收函数
    pArr = zeros(nsta, 1);     % 存储平均射线参数
    stlas = zeros(nsta, 1);    % 台站纬度
    stlos = zeros(nsta, 1);    % 台站经度
end

%% 子函数：加载速度模型
function [z, r, vp, vs] = loadVelocityModel()
    % 加载 AK135 速度模型
    [z, r, vp, vs, ~, ~] = ak135('cont');
end

%% 子函数：叠加波形
function [seis, pArr] = stackWaveforms(gather, dz, zmax, z, vp, vs)
    % 检查所有记录长度是否一致
    rfsAll = cellfun(@(rf) rf.itr, {gather.RF}, 'UniformOutput', false);
    timeAll = cellfun(@(rf) rf.ittime, {gather.RF}, 'UniformOutput', false);

    % 提取射线参数
    raypAll = cellfun(@(travelinfo) travelinfo.rayParam / 6371, {gather.TravelInfo}, 'UniformOutput', false);
    raypAll = cell2mat(raypAll);

    % 时间-深度转换
    [~, rfsAll_depth, ~] = rf_migrate(timeAll, rfsAll, raypAll, dz, zmax, z, vp, vs);

    % 叠加波形
    seis = mean(cell2mat(rfsAll_depth), 2, 'omitnan');
    pArr = mean(raypAll, 'omitnan');
end

%% 子函数：平滑地震数据
function seisout = smoothSeismicData(seisCell,ngrid_x,ngrid_y)
    % 将 cell 转换为矩阵
    seisout = cell2mat(seisCell);

    % 平滑处理
    kernel = ones(ngrid_x, ngrid_y) / (ngrid_x * ngrid_y);
    seisout = conv2(seisout, kernel, 'same');
end

%% 子函数：绘制叠加剖面
function plotStackedSection(seisout, depth0, dmoho)
    figure('Name', 'Stacked RF section', 'Color', 'w', 'Position',[10 10 800 400]);
    imagesc(1:size(seisout, 2), depth0, seisout); hold on;
    colormap(seismic(1)); caxis([-0.05 0.05]);
    plot(1:size(seisout, 2),dmoho,'r--'); hold off;
    ylim([0 100]);
    xlabel('Station index'); ylabel('Depth (km)');
    colorbar;
    set(gca,'fontsize',14)
end

%% 子函数：提取莫霍面深度
function [dmoho, amoho] = extractMohoDepth(seisout, depth0)
    drange = [40, 60];  % 莫霍面深度范围
    depthIndex = find(depth0 >= drange(1) & depth0 <= drange(2));
    nx = size(seisout, 2);
    dmoho = zeros(nx, 1);
    amoho = zeros(nx, 1);

    for i = 1:nx
        ampData = seisout(depthIndex, i);
        rmsVal = rms(ampData);
        [peakVal, peakInd] = max(ampData);

        if peakVal < 1.0 * rmsVal
            dmoho(i) = NaN;  % 峰值不明显，标记为 NaN
        else
            indMoho = depthIndex(peakInd);
            dmoho(i) = depth0(indMoho);
            amoho(i) = peakVal;
        end
    end
end

%% 子函数：剔除离群值并插值
function [Fmoho, idx] = interpolateMohoDepth(stlas, stlos, dmoho)
    % 剔除离群值
    mu = mean(dmoho, 'omitnan');
    sig = std(dmoho, 'omitnan');
    idx = (dmoho >= mu - 3 * sig) & (dmoho <= mu + 3 * sig);

    % 散点插值
    Fmoho = scatteredInterpolant(stlos(idx), stlas(idx), dmoho(idx), 'natural', 'none');
end

%% 子函数：绘制莫霍面深度图
function plotMohoMap(stlas, stlos, dmoho, Fmoho)
    % plotMohoMap 使用 MATLAB Mapping Toolbox 绘制 Moho 深度等值线图
    % stlas, stlos: 台站（或观测点）的纬度和经度
    % dmoho: 原始传入的 Moho 深度向量（如果用不到可以忽略）
    % Fmoho: Moho 深度插值函数，比如 TriScatteredInterp 或者 scatteredInterpolant

    figure('Name', 'Moho depth map', 'Color', 'w', 'Position',[10 10 800 400]);

    % 1. 计算绘图范围
    latlim = [min(stlas), max(stlas)];
    lonlim = [min(stlos), max(stlos)];

    % 2. 构建规则网格，用于绘制等值线
    [lonGrid, latGrid] = meshgrid(lonlim(1):0.05:lonlim(2), ...
                                  latlim(1):0.05:latlim(2));
    % 从插值函数获得此网格上的 Moho 值
    VI = Fmoho(lonGrid, latGrid);

    % 3. 设置地图投影，并指定经纬度范围
    %    axesm 中常见 'MapProjection' 选项有 'lambertstd', 'lambert', 'lambertconic'
    %    也可以考虑用 worldmap(latlim, lonlim) 来简化。
    ax = worldmap(latlim, lonlim);
%     ax = axesm('lambertstd', ...
%                'MapLatLimit', latlim, ...
%                'MapLonLimit', lonlim, ...
%                'FLineWidth', 2, ...    % 地图框线宽
%                'FontSize', 12);       % 地图文字大小
    hold on;

    % 4. 绘制等值线 (contourfm/contourm)
    %    注意 contourfm 的输入顺序是 (lat, lon, Z)。因为 meshgrid 我们先生成 (lonGrid, latGrid)，
    %    所以这里要用 contourfm(latGrid, lonGrid, VI)。
    contourfm(latGrid, lonGrid, VI, 15);  % 绘制 10 条等值线
    % 调整色表：先用 jet，再 flipud 翻转
    colormap(flipud(jet));
    colorbar('Location', 'eastoutside', 'FontSize', 12);

    % 5. 在地图上绘制测点散点
    plotm(stlas, stlos, 'k^', 'MarkerSize', 6, 'LineWidth', 1);

    % 6. 绘制海岸线/行政区线
    %    这里举例用 landareas.shp。你可以用其它矢量数据替代。
    %    如果你有 gshhs_i.shp，也可以改为 geoshow('gshhs_i.shp', 'Color','k')
    try
        geoshow('landareas.shp', 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
    catch
        % 如果没有 landareas.shp 文件，可以用内置 coast.mat 的数据:
        % load coast
        % plotm(lat, long, 'k', 'LineWidth',1);
        warning('找不到 landareas.shp，使用内置海岸线数据替代。');
        load coast
        plotm(lat, long, 'k', 'LineWidth',1);
    end

    % 7. 打开经纬度网格
    setm(ax, 'FontSize', 14, ...
             'Grid', 'on', ...         % 显示网格
             'Frame', 'on', ...        % 显示地图边框
             'MeridianLabel', 'on', ...% 标注经度
             'ParallelLabel', 'on', ...% 标注纬度
             'MLineLocation', 0.5, ...   % 经度网格线间隔
             'PLineLocation', 0.5);      % 纬度网格线间隔

    mlabel('on');  % 打开经度标注
    plabel('on');  % 打开纬度标注    
    gridm('GLineStyle', '--', 'GColor', 'k', 'GLineWidth', 0.5);  % 设置网格线样式

    title('Moho depth', 'FontSize', 16);

    % 8. 绘制 Baiyan Ebo 矿区矩形
    latB = [41.65, 41.65, 41.8833, 41.8833, 41.65];
    lonB = [109.7833, 110.0667, 110.0667, 109.7833, 109.7833];
    plotm(latB, lonB, 'b', 'LineWidth', 1);

    % 9. 读取并绘制断层
    fault_files = dir('./visualization/faults/*txt');
    for k = 1:length(fault_files)
        faults = read_faults_gmt(fullfile(fault_files(k).folder, fault_files(k).name));
        for l = 1:length(faults)
            fault = faults{l};
            % fault(:,1) -> lon, fault(:,2) -> lat
            geoshow(fault(:,2), fault(:,1),'DisplayType','line','LineWidth',2,'Color','k');
        end
    end

end

%% 子函数：创建莫霍面结构体
function mohoStruct = createMohoStruct(stlas, stlos, dmoho, amoho)
    mohoStruct.lat = stlas;
    mohoStruct.lon = stlos;
    mohoStruct.dmoho = dmoho;
    mohoStruct.amoho = amoho;
end