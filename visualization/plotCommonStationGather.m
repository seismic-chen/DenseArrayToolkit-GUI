function plotCommonStationGather(DataStruct, station)
    % plotCommonStationGather  绘制同一台站的RF道集
    %
    % 参数说明：
    %   DataStruct : 数据结构体数组，每个元素包含 TravelInfo、RF、StationInfo 等
    %   station    : 字符串，指定要绘制的台站名

    %========== 寻找所有满足 station 条件的索引 ==========
    stationList = cellfun(@(stationinfo) stationinfo.sta, {DataStruct.StationInfo}, 'UniformOutput', false);
    idx = find(strcmp(stationList, station));
    
    % 若没有找到对应的台站数据，可以给出警告/提示
    if isempty(idx)
        warning('未找到名为 %s 的台站数据！', station);
        return;
    end

    %========== 统一可配置参数 ==========
    figName        = ['Common Station Gather: ', station];
    figPosition    = [500 500 600 1200];
    scaling_factor = 10;       % 纵向放大倍数
    thr            = 0.01;     % 阈值
    timeRange      = [0, 30];  % 横坐标时间范围
    distRange      = [25,100]; % 纵坐标距离范围
    
    %========== 创建图形与基本设置 ==========
    figure('Name',figName,'Position',figPosition,'Color','w');
    hold on;  % 只需在此处开启一次
    title(figName, 'Interpreter', 'none');  % 避免特殊字符导致的问题
    
    %========== 绘制每条道的RF信息 ==========
    for n = idx
        % 取出该道需要的数据
        x = DataStruct(n).RF.itr * scaling_factor;
        t = DataStruct(n).RF.ittime;
        dist = DataStruct(n).TravelInfo.distDeg;

        % 1) 绘制原始波形
        p = plot(t, x + dist, 'k', 'LineWidth', 0.7);

        % 2) 正极性部分着色
        fillPositive(t, x, dist, thr); hold on;

        % 3) 负极性部分着色
        fillNegative(t, x, dist, thr); hold on;
    end
    hold off;

    %========== 设置坐标轴与标签 ==========
    xlim(timeRange);
    ylim(distRange);
    xlabel('Time (sec)');
    ylabel('Distance (deg)');
    set(gca, 'FontSize', 14);
end

%===================== 辅助函数：正极性填充 =====================%
function fillPositive(t, x, dist, thr)
    % fillPositive 对 x>thr 的部分进行填充着色
    upper               = x;
    upper(upper <= thr) = 0;
    upper(1)            = 0;
    upper(end)          = 0;
    
    lower = zeros(size(x));
    
    % jbfill 的参数含义(根据版本可能略有差异)：
    % jbfill(x轴, 上边界, 下边界, 面颜色, 边框颜色, 填充方式，透明度)
    jbfill(t, upper + dist, lower + dist, 'r', 'none', 1, 0.5);
end

%===================== 辅助函数：负极性填充 =====================%
function fillNegative(t, x, dist, thr)
    % fillNegative 对 x<-thr 的部分进行填充着色
    lower               = x;
    lower(lower >= -thr)= 0;
    lower(1)            = 0;
    lower(end)          = 0;
    
    upper = zeros(size(x));
    
    jbfill(t, upper + dist, lower + dist, [0.17,0.17,0.17], 'none', 1, 0.5);
end
