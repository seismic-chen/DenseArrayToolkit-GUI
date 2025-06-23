function plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profile)
%% create 2D cross-section
% 创建插值函数
F = scatteredInterpolant(X(:),Y(:),Z(:),V(:));
depth0 = 0:0.5:max(Z(:)); % 定义深度范围

% 加载颜色映射表
cmap = load('roma.mat');
% 加载DEM数据
dem = load('Qaidam_DEM.mat');
demX = X(:,:,1);
demY = Y(:,:,1);
[nx,ny] = size(demX);
[LON, LAT] = ProjectedCoordsTolatlon(demX(:),demY(:), gridStruct);
demZ = interp2(dem.demLon,dem.demLat,dem.demZ,LON,LAT)/1000.;
demZ = reshape(demZ,nx,ny);

% 创建图形窗口
figure;
set(gcf,'Position',[0 0 1000 1000],'Color','w')
ax1 = subplot(1,1,1);
hold on;
% 绘制三维切片图
h = slice(X,Y,Z,V,[],[],100);
colormap(ax1, flipud(cmap.roma));  % 为波形数据设置不同的色标
cmax = rms(V(:));  % 计算色标的最大值
caxis([-cmax, cmax]);  % 设置波形数据的色标范围

xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
set(gca,'FontSize',16)

% 初始化变量
VAll=[];
distAll = [];
depthAll = [];
totalDist = 0;
elevAll = [];

% 遍历每个剖面
for n = 1:length(profile)
    x1 = profile{n}(1,1);
    y1 = profile{n}(1,2);
    x2 = profile{n}(2,1);
    y2 = profile{n}(2,2);
    l = (y2-y1)/(x2-x1);
    npoints = 100;
    xp = linspace(x1,x2,npoints)';
    yp = y1+(xp-x1)*l;
    dist0 = sqrt((xp-x1).^2 + (yp-y1).^2);
    % 绘制RFs的几何形状
    Y_profile = repmat(yp',length(depth0),1);
    X_profile = repmat(xp',length(depth0),1);
    Z_profile = repmat(depth0',1,length(xp));
    DIST_profile = repmat(dist0',length(depth0),1);
    Vprofile = F(X_profile,Y_profile,Z_profile);
    Vprofile(isnan(Vprofile)) = NaN;

    scatter(xp,yp,50,'ko','filled');
    surface(X_profile,Y_profile,Z_profile,Vprofile,'EdgeColor','none')

    colormap(flipud(cmap.roma))
    cmax = rms(Vprofile(:));
    caxis([-cmax,cmax]);
    
    grid on; box on;
    view(140,20)
    if n > 1
        totalDist = max(distAll(:));
    end
    distAll = [distAll, DIST_profile+totalDist];
    depthAll = [depthAll,Z_profile];
    VAll = [VAll, Vprofile];
    % 提取DEM文件中的地表高程
    [lonp, latp] = ProjectedCoordsTolatlon(xp, yp, gridStruct);
    elevp = interp2(dem.demLon,dem.demLat,dem.demZ,lonp,latp)/1000.;
    elevAll = [elevAll; elevp];
end
zlim([-30 100])
distElev = distAll(1,:)';

% 创建第二个坐标轴，用于绘制DEM
ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...       % 背景透明
           'XTick', [], 'YTick', [], 'ZTick', [], ...
           'View', ax1.View); % 隐藏坐标刻度
% linkaxes([ax1,ax2])
linkprop([ax1,ax2],{"CameraPosition","CameraUpVector"});
rotate3d on
hdem = surface(ax2,demX,demY,-5*demZ);
colormap(ax2,'parula')
set(hdem,'EdgeColor','none','FaceAlpha',0.8)
zlim([ax1.ZLim])
xlim([ax1.XLim])
ylim([ax1.YLim])
set(ax2,'ZDir','reverse','Visible','off')

% % 创建新的图形窗口，用于绘制剖面图
% fig = figure();
% set(gcf,'Position',[50 50 1200 600],'Color','w')
% subplot(511,'Parent',fig)
% plot(distElev,elevAll,'k','linewidth',2);
% patch([distElev(1);distElev;distElev(end)],[0;elevAll;0],[0.8 0.8 0.8]);
% ylabel('Elev. (km)')
% xlim([0 max(distAll(1,:))]);
% set(gca,'fontsize',14)
% ylim([1 6]);
% subplot(5,1,2:5,'Parent',fig);
% h = pcolor(distAll,depthAll,VAll);
% ylim([0 100]);
% colormap(flipud(cmap.roma))
% cmax = rms(VAll(:));
% caxis([-cmax,cmax]);
% set(h,'EdgeColor','none')
% axis ij
% xlabel('Distance (km)')
% ylabel('Depth (km)')
% set(gca,'fontsize',14)