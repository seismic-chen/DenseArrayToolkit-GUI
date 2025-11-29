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
    load('./visualization/colormap/roma.mat')
    station = getStations(DataStruct); 
    stationList = {station.sta};
    stla = [station.stla]; stlo = [station.stlo];
    nsta = length(stationList);
    arraytype = station(1).network;
    if strcmp(arraytype , 'SL')==1
        load('./visualization/SichuanLongmenshan.mat');
    elseif strcmp(arraytype , 'BY') == 1
        load('./visualization/Baiyanebo_DEM.mat');
    end

    F = scatteredInterpolant(demLon(:), demLat(:), demZ(:), 'nearest');
    staelev = F(stlo, stla)/1000;

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
    plotMohoMap(arraytype,stlas, stlos, dmoho, Fmoho);

    % 10) 绘制叠加剖面
    plotStackedSection(arraytype,seisout, depth0, dmoho, staelev,roma);

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
function plotStackedSection(arraytype,seisout, depth0, dmoho, staelev,roma)
    if arraytype == 'SL'
        seisout = fliplr(seisout);
        dmoho = fliplr(dmoho');
        staelev = fliplr(staelev);

        
    elseif arraytype == 'BY'

    else

    end
    f1 = figure('Name', 'Stacked RF section', 'Color', 'w', 'Position',[10 10 800 450]);
    axRF = axes('Parent',f1,'Position', [0.1 0.12 0.75 0.6]);
    rf_pos = axRF.Position;
    imagesc(1:size(seisout, 2), depth0, seisout); hold on;
    colormap(flipud(roma)); clim([-0.05 0.05]);
    plot(1:size(seisout, 2), dmoho, 'r--', 'LineWidth', 0.8); hold off;
    ylim([0 100]);xlim([0 size(seisout, 2)])
    ylabel('Depth (km)');xlabel('# of station')
    set(axRF, 'FontSize', 18,'LineWidth',1.5)
    title(axRF, 'Stacked RFs');

    c = colorbar(axRF,'Location','eastoutside');
    c.Position(1) = rf_pos(1) + rf_pos(3) + 0.01;
    c.Position(2) = rf_pos(2);
    c.Position(3) = 0.02;   c.Position(4) = rf_pos(4);

    axRF.Position = rf_pos;
    gap   = 0.08;  hTop  = 0.18;
    yTop  = rf_pos(2) + rf_pos(4) + gap;

    axElev = axes('Parent',f1,'Position', [rf_pos(1) yTop rf_pos(3) hTop]);
    plot(1:size(seisout, 2), staelev, 'b-', 'LineWidth', 2);
    xlim([0 size(seisout, 2)]);ylim([0 5])
    grid on; box on;
    xlabel('Station index'); ylabel('Elev. (km)');
    axElev.XTickLabel = [];   axElev.XLabel.String = '';
    set(axElev,'FontSize',18,'LineWidth',1.5);
    linkaxes([axRF, axElev], 'x');
end

function [dmoho, amoho] = extractMohoDepth(seisout, depth0)
    drange = [40, 60];  
    depthIndex = find(depth0 >= drange(1) & depth0 <= drange(2));
    nx = size(seisout, 2);
    dmoho = zeros(nx, 1);
    amoho = zeros(nx, 1);

    for i = 1:nx
        ampData = seisout(depthIndex, i);
        rmsVal = rms(ampData);
        [peakVal, peakInd] = max(ampData);

        if peakVal < 1.0 * rmsVal
            dmoho(i) = NaN; 
        else
            indMoho = depthIndex(peakInd);
            dmoho(i) = depth0(indMoho);
            amoho(i) = peakVal;
        end
    end
end

function [Fmoho, idx] = interpolateMohoDepth(stlas, stlos, dmoho)
    mu = mean(dmoho, 'omitnan');
    sig = std(dmoho, 'omitnan');
    idx = (dmoho >= mu - 3 * sig) & (dmoho <= mu + 3 * sig);

    Fmoho = scatteredInterpolant(stlos(idx), stlas(idx), dmoho(idx), 'natural', 'none');
end

function plotMohoMap(arraytype,stlas, stlos, dmoho, Fmoho)

    latlim = [min(stlas), max(stlas)];
    lonlim = [min(stlos), max(stlos)];
    [lonGrid, latGrid] = meshgrid(lonlim(1):0.05:lonlim(2), latlim(1):0.05:latlim(2));
    VI = Fmoho(lonGrid, latGrid);

    figure('Name', 'Moho depth map', 'Color', 'w', 'Position',[10 10 800 400]);
    contourfm(latGrid, lonGrid, VI, 15);
    box on;colormap(flipud(turbo));
    title('Moho depth')
    xlabel('Longitude (deg)');ylabel('Latitude (deg)');
    set(gca,'FontSize', 18,'BoxStyle','full','LineWidth',1.5);
    hold on
    plot( stlos,stlas, 'k^', 'MarkerSize', 6, 'LineWidth', 1);
    axis tight;
    colorbar
    latB = [41.65, 41.65, 41.8833, 41.8833, 41.65];
    lonB = [109.7833, 110.0667, 110.0667, 109.7833, 109.7833];
    plot(lonB, latB,'k', 'LineWidth', 0.8,'linestyle','--')
    
    if arraytype == 'BY'
        fault_files = dir('./visualization/faults/*txt');
        for k = 1:length(fault_files)
            faults = read_faults_gmt(fullfile(fault_files(k).folder, fault_files(k).name));
            for l = 1:length(faults)
                fault = faults{l};
                % fault(:,1) -> lon, fault(:,2) -> lat
                geoshow(fault(:,2), fault(:,1),'DisplayType','line','LineWidth',2,'Color','k');
            end
        end
        axis([lonlim,latlim,])
    else

    end
    hold off;

end

function mohoStruct = createMohoStruct(stlas, stlos, dmoho, amoho)
    mohoStruct.lat = stlas;
    mohoStruct.lon = stlos;
    mohoStruct.dmoho = dmoho;
    mohoStruct.amoho = amoho;
end