function plotWaveforms(DataStruct, trace_index)
% plotWaveforms  绘制三分量地震记录与接收函数

%========= 获取数据 =========%
t       = DataStruct(trace_index).TimeAxis.t_resample;
T       = DataStruct(trace_index).Waveforms.dataProcessed(:,1);
R       = DataStruct(trace_index).Waveforms.dataProcessed(:,2);
Z       = DataStruct(trace_index).Waveforms.dataProcessed(:,3);
ittime  = DataStruct(trace_index).RF.ittime;
itr     = DataStruct(trace_index).RF.itr;
wltime  = DataStruct(trace_index).RF.wltime;
wlr     = DataStruct(trace_index).RF.wlr;

%========= 统一图像设置 =========%
figure('Position',[10,10,1000,800],'Color','w');

% 如果想对前三个波形统一给出显示区间，可以这样定义：
timeRange3C = [t(1), t(1) + 120];  % 示例：显示前120s

%========= 绘制 3 分量 =========%
% 1) R 分量
subplot(5,1,1);
plot(t, R, 'k');
xlim(timeRange3C);
ylabel('R');
title('Three Components Seismogram');

% 2) T 分量
subplot(5,1,2);
plot(t, T, 'k');
xlim(timeRange3C);
ylabel('T');

% 3) Z 分量
subplot(5,1,3);
plot(t, Z, 'k');
xlim(timeRange3C);
ylabel('Z');

%========= 绘制接收函数 =========%
% 设置统一阈值
thr = 0.01;

% 4) wlr
subplot(5,1,4);
thresholdFillPlot(wltime, wlr, thr);
xlim([wltime(1), 30]);
ylim([-0.25, 0.5]);
ylabel('RF (wlr)');

% 5) itr
subplot(5,1,5);
thresholdFillPlot(ittime, itr, thr);
xlim([ittime(1), 30]);
ylim([-0.25, 0.5]);
xlabel('Time (sec)');
ylabel('RF (itr)');

end

function thresholdFillPlot(t, x, thr)
% thresholdFillPlot 在同一张图上画原始曲线，并对超过阈值的部分用颜色填充

plot(t, x, 'k');
hold on;

% 将高于阈值的部分放入 upper，其余置0
upper          = x;
upper(upper<=thr) = 0;
% 为保证填充区域闭合，首尾也可置零
upper(1)       = 0;
upper(end)     = 0;

lower = zeros(size(x));
% 这里使用 jbfill，如果你使用的是其他填充函数，请根据需要替换
jbfill(t, upper, lower, 'r', 'none', 1, 1.0);

hold off;
end
