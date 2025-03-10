 % Sep. 4, 2017, Yunfeng Chen, Global Seismology Group, University of
% Alberta
% convert x y to lat lon
% function [lon,lat] = xy2latlon(x,y,lonmin,latmin)
function [lon,lat] = xy2latlon(x,y,lonmin,latmin)

dconst=6378.0*pi/180.0;
lat=y/dconst+latmin;
lon=x./dconst./cos(lat/180.*pi)+lonmin;