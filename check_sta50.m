
id = '20250515005052';
rfmatrix_idx = zeros(1500,238);
Rcomp_idx0 = rfmatrix_idx;
Zcomp_idx0 = rfmatrix_idx;


for idx = foundIndices
    stanm = gather(idx).StationInfo.sta;
    tmp = str2double(strsplit(stanm,stanm(1:2)));
    stanum = tmp(2);
    rfmatrix_idx(:, stanum) = gather(idx).RF.itr;

    Rcomp_idx0(:,stanum) = gather(idx).Waveforms.dataProcessed(:,2); %  R
    Zcomp_idx0(:,stanum) = gather(idx).Waveforms.dataProcessed(:,3); %  Z

end

%% 

dt = 0.1;
Zcomp_idx = bandpassSeis(Zcomp_idx0, dt, 0.2, 2, 3);
Rcomp_idx = bandpassSeis(Rcomp_idx0, dt, 0.2, 2, 3);

% Zcomp_idx = Zcomp_idx(250:480,:);
% Rcomp_idx = Rcomp_idx(250:480,:);
% 
% Zcomp_idx(1:295,:)=eps;
% Rcomp_idx(1:295,:)=eps;

% Zcomp_idx(700:end,:)=eps;
% Rcomp_idx(700:end,:)=eps;

nt=size(Zcomp_idx,1);
phaseshift=5;
tt= (0:nt-1)*dt;

gauss = 2.5;
itmax = 100;
minderr = 1e-5;

for ii = 1:238
    [rfs, itrms] = makeRFitdecon_la_norm( ...
        Rcomp_idx(:,ii), Zcomp_idx(:,ii), dt, nt, phaseshift, gauss, itmax, minderr);

    rfmatrix_idx(1:nt,ii) = rfs;
end

rf50 = rfmatrix_idx(:,50);
figure('Position',[200 200 800 300],'Color','w');
hold on;
plot(trf,rf50,'LineWidth',1.5,'Color','r');
xlim([-2 20]);
xlabel('Time (sec)');box on;
set(gca,'fontsize',18,'linewidth',1.2)

n1 = 40;n2 = 60;
figure('Position',[200 200 1200 320],'Color','w')
subplot(131)

for jj = n1:n2
    if jj == 50 || jj == 46 || jj == 56 || jj == 57
        plot(jj+Zcomp_idx(:,jj)*50,tt,'Color','r','LineWidth',0.8)
    elseif jj == 44
        continue
    else
        plot(jj+Zcomp_idx(:,jj)*50,tt,'Color',[0.3 0.3 0.3],'LineWidth',0.8)
    end
    
    hold on
end
title('Z component')
xlabel('Station');ylabel('Time (s)')
% ylim([25 45])
box on
set(gca,'fontsize',18,'YDir','reverse')

subplot(132)
for jj = n1:n2
    if jj == 50 || jj == 46 || jj == 56 || jj == 57
        plot(jj+Rcomp_idx(:,jj)*50,tt,'Color','r','LineWidth',0.8)
    elseif jj == 44
        continue
    else
        plot(jj+Rcomp_idx(:,jj)*50,tt,'Color',[0.3 0.3 0.3],'LineWidth',0.8)
    end
    
    hold on
end
title('R component')
% ylim([25 45])
xlabel('Station');ylabel('Time (s)')
box on
set(gca,'fontsize',18,'YDir','reverse')

subplot(133)
wigb(rfmatrix_idx(:,n1:n2),2,n1:n2,trf,[],[7,11,17,18]);
title('RF')
ylim([-2 20])
box on
xlabel('Station');ylabel('Time (s)')
set(gca,'fontsize',18)






% 
% 
% 
% 
% 
% figure('Position',[200 200 1200 320],'Color','w')
% subplot(131)
% wigb(Zcomp(:,39:59),1.2,distArr(39:59),trf+5,[],11);
% title('Z component')
% xlabel('Distance (deg)');ylabel('Time (s)')
% ylim([25 45])
% box on
% set(gca,'fontsize',18)
% subplot(132)
% wigb(Rcomp(:,39:59),1.2,distArr(39:59),tt,[],11);
% title('R component')
% xlabel('Distance (deg)');ylabel('Time (s)')
% ylim([25 45])
% box on
% set(gca,'fontsize',18)
% subplot(133)
% wigb(rfmatrix(:,39:59),2,distArr(39:59),trf,[],11);
% title('RF')
% ylim([-2 20])
% box on
% xlabel('Distance (deg)');ylabel('Time (s)')
% set(gca,'fontsize',18)
