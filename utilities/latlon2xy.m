% Sep. 4, 2017, Yunfeng Chen, Global Seismology Group, University of
% Alberta
% convert lat lon to x y
function [x,y] = latlon2xy(lon,lat,lonmin,latmin)
dconst=6378.0*pi/180.0;
y = (lat-latmin)*dconst;
x = (lon-lonmin).*cos(lat/180*pi)*dconst;                              