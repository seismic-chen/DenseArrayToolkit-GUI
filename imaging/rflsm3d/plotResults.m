function plotResults(mig,pre_rfm,lsmig,pre_rflsm,param,mask)
    load roma;
    param = param.paramMig;
% paramters
    x = param.x ;
    y = param.y ;
    z = param.z;

%% plot migration images
    yc = (max(y)+min(y))/2;
    xc = (max(x)+min(x))/2;
    
    nnx = find(y < yc + 6 & y > yc - 6,1);    % y axis  km
    nny = find(x < xc + 5 & x > xc - 5,1);    % x axis  km

    if isempty(nnx) || isempty(nny)
        nnx = round(numel(y)/2);
        nny = round(numel(x)/2);
    end

    stax = x(mask(:,nnx)==1);
    stay = y(mask(nny,:)==1);
    
    figure
    set(gcf,'Position',[1,200,1950,700],'Color','white')
    subplot(2,5,[1 2 3])
    imagesc(x,z,mig(:,:,nnx)./max(mig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('Migration image')
    caxis([-0.1 0.1]);
    axis equal
    colormap(flipud(roma))
    ylim([0 100])
    xlim([min(x) max(x)])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('X (km)')
    ylabel('Depth (km)')
    hold on
    scatter(stax,zeros(1,numel(stax)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    
    hold off
    text(-15,-12,'(a)','FontSize',20)
    
    
    subplot(2,5,[6 7 8])
    imagesc(x,z,lsmig(:,:,nnx)./max(lsmig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('LSM image')
    axis equal
    colormap(flipud(roma))
    ylim([0 100])
    xlim([min(x) max(x)])
    colormap(roma)
    caxis([-0.1 0.1])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('X (km)')
    ylabel('Depth (km)')
    hold on
    scatter(stax,zeros(1,numel(stax)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);

    hold off
    text(-15,-12,'(c)','FontSize',20)
    
    subplot(2,5,[4 5])
    imagesc(y,z,squeeze(mig(:,nny,:))./max(mig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('Migration image')
    caxis([-0.1 0.1]);
    axis equal
    colormap(flipud(roma))
    ylim([0 100])
    xlim([min(y) max(y)])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('Y (km)')
    ylabel('Depth (km)')
    hold on
    scatter(stay,zeros(1,numel(stay)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    hold off
    text(-15,-12,'(b)','FontSize',20)
    
    subplot(2,5,[9 10])
    imagesc(y,z,squeeze(lsmig(:,nny,:))./max(lsmig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('LSM image')
    axis equal
    colormap(flipud(roma))
    ylim([0 100])
    xlim([min(y) max(y)])
    caxis([-0.1 0.1])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('Y (km)')
    ylabel('Depth (km)')
    hold on
    scatter(stay,zeros(1,numel(stay)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    hold off
    text(-15,-12,'(d)','FontSize',20)


end