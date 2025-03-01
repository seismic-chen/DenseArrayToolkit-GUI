function gather = project_stations(gather,center,direction,if_plot)
    % PROJECT_STATIONS 将台站位置投影到测线
    %
    % 输入参数:
    %   gather - Struct of common event gather
    %   center - center of the best fit profile
    %   direction - direction of the best fit profile
    %   if_plot = if plot figure or not
    % 输出参数:
    %   gather - Struct of common event gather
    %
    % 绘制图形:
    %   1. 台站位置
    %   2. 最佳拟合测线
    %   3. 台站位置在测线上的投影点
    
    if nargin < 3
        if_plot = 0;
    end
    % 将经度和纬度合并为一个矩阵
    slat = cellfun(@(stationinfo) stationinfo.stla, {gather.StationInfo}, 'UniformOutput', false);
    slon = cellfun(@(stationinfo) stationinfo.stlo, {gather.StationInfo}, 'UniformOutput', false);
    slat = cell2mat(slat)';
    slon = cell2mat(slon)';
    data = [slon(:), slat(:)];
    centered_data = data - center;
    projected_points = center + (centered_data * direction) * direction';
    slatp=projected_points(:,2);
    slonp=projected_points(:,1);

    % 生成测线上的点
    t = linspace(-2, 2, 100); % 根据需要调整范围
    line_points = center + t' * direction';
    lat1 = line_points(1,2);
    lon1 = line_points(1,1);
    % calculate the receiver location along the profile
    [deg0,~]= distance(lat1,lon1,slatp,slonp);
    rx = deg0*2*pi*6371/360;

    for n=1:length(gather)
        gather(n).RF.rx = rx(n);
    end
    if if_plot
        % 绘制原始台站位置
        figure;
        plot(data(:, 1), data(:, 2), 'bo', 'DisplayName', '台站位置');
        hold on;

        % 绘制最佳拟合测线
        plot(line_points(:, 1), line_points(:, 2), 'r-', 'DisplayName', '最佳拟合测线');

        % 绘制投影点
        plot(projected_points(:, 1), projected_points(:, 2), 'gx', 'DisplayName', '投影点');
        
        % 连接原始点和投影点
        for i = 1:size(data, 1)
            plot([data(i, 1), projected_points(i, 1)], [data(i, 2), projected_points(i, 2)], 'k--');
        end

        % 设置图形属性
        xlabel('经度');
        ylabel('纬度');
        title('台站位置及其在最佳拟合测线上的投影');
%         legend show;
        axis equal;
        grid on;
        hold off;
    end
end