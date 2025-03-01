function plotMigrationResults(dmig,dmigls,x,z)
% dmig = migResult.mig;
% dmigls = migResult.migls;
% x = migResult.x;
% z = migResult.z;

d2d=mean(dmig,3);
d2dls=mean(dmigls,3);
xmax = max(x);
figure();
set(gcf,'Position',[100 100 800 800],'color','w')
subplot(211)
imagesc(x,z, d2d/max(d2d(:))); hold on;

axis([0 xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('Migration image')
set(gca,'fontsize',14)
cmax=0.3;
caxis([-cmax cmax]);
colorbar
text(-0.12,0.98,'a)','Units','normalized','FontSize',18)

subplot(212)
imagesc(x,z, d2dls/max(d2dls(:))); hold on;
axis([0 xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('LSM')
set(gca,'fontsize',14)
caxis([-cmax cmax]);
colorbar
colormap(seismic(1));
text(-0.12,0.98,'b)','Units','normalized','FontSize',18)
