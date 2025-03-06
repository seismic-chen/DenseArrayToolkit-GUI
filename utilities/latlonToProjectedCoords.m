function [rx, ry] = latlonToProjectedCoords(stlo, stla, gridStruct)
    % 将经度和纬度合并为一个矩阵
%     stla = cellfun(@(stationinfo) stationinfo.stla, {gather.StationInfo}, 'UniformOutput', false);
%     stlo = cellfun(@(stationinfo) stationinfo.stlo, {gather.StationInfo}, 'UniformOutput', false);
%     stlo = cell2mat(stlo)';
%     stla = cell2mat(stla)';

    % 将经纬度转换为笛卡尔坐标系
    [stationX, stationY] = latlon2xy(stlo, stla, gridStruct.originLon, gridStruct.originLat);
    
    % 执行PCA分析
    coords = [stationX(:), stationY(:)];
    projection_on_principal_axis = coords * gridStruct.coeff(:, 1);  % 主轴方向投影
    projection_on_secondary_axis = coords * gridStruct.coeff(:, 2);  % 次轴方向投影

    % 台站投影位置
    rx = projection_on_principal_axis; % 主轴上的投影
    ry = projection_on_secondary_axis; % 次轴上的投影

end