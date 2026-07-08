#!/bin/bash
# GMT 6.4.0 兼容版本
# 定义区域
region=105/115/38/45
grd_file=earth_relief_15s.grd

# 下载高分辨率的地表高程数据
# 添加 -r 参数指定gridline registration
gmt grdcut @earth_relief_15s -R$region -G$grd_file -r

# 如果需要查看数据范围
gmt grdinfo $grd_file
