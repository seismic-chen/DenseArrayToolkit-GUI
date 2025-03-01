function [line_points, projected_points, mean_data, direction] = fit_and_project_stations(slon,slat,profile_length,if_plot)
    % FIT_AND_PROJECT_STATIONS 拟合台站位置到最佳测线并投影
    %
    % 输入参数:
    %   slon - 台站的经度 (列向量)
    %   slat - 台站的纬度 (列向量)
    %   profile_length - the length of profile in degree
    %   if_plot = if plot figure or not
    % 输出参数:
    %   line_points - 最佳拟合测线上的点 (用于绘图)
    %   projected_points - 台站位置在测线上的投影点
    %
    % 绘制图形:
    %   1. 台站位置
    %   2. 最佳拟合测线
    %   3. 台站位置在测线上的投影点
    
    if nargin < 3
        if_plot = 0;
    end
    % 将经度和纬度合并为一个矩阵
    data = [slon(:), slat(:)];

    % 对数据进行中心化
    mean_data = mean(data);
    centered_data = data - mean_data;

    % 使用PCA找到最佳拟合测线
    [coeff, ~, ~] = pca(centered_data);

    % 提取第一个主成分方向（最佳拟合测线的方向）
    direction = coeff(:, 1);

    % 生成测线上的点
    t = linspace(-profile_length/2, profile_length/2, 100); % 根据需要调整范围
    line_points = mean_data + t' * direction';

    % 将台站位置投影到测线上
    projected_points = mean_data + (centered_data * direction) * direction';
    
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