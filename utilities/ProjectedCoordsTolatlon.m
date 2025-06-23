function [lon, lat] = ProjectedCoordsTolatlon(rx, ry, gridStruct)
OriginalCoord = [rx, ry]* gridStruct.coeff'; 
rxInOriginalCoord = OriginalCoord(:,1);

ryInOriginalCoord = OriginalCoord(:,2);
% 将笛卡尔坐标系转换为经纬度
[lon, lat] = xy2latlon(rxInOriginalCoord, ryInOriginalCoord, gridStruct.originLon, gridStruct.originLat);