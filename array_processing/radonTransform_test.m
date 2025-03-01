function radonTransform(DataStruct)
lows  = 0.1;    % Low corner frequency (Hz)
highs = 1.2;    % High corner frequency (Hz)

%% Prepare common shot gather geometry
eventlist = getEvent(DataStruct);
%% Process events
for i = 1:length(eventlist)
    CommonEventGather = getCommonEventGather(DataStruct,eventlist{i});
    % Criteria: at least 60 traces
    if length(CommonEventGather) < 60
        continue
    end

    % Display progress
    disp(['Processing event# ', num2str(i)]);

    d_z = [];
    d_r = [];
    for n = 1:length(CommonEventGather)
        dataProcessed = CommonEventGather(n).Waveforms.dataProcessed;
        d_z = [d_z, dataProcessed(:,3)];
        d_r = [d_r, dataProcessed(:,2)];
    end

    t        = CommonEventGather(1).RF.ittime;
    dt       = t(2)-t(1);

    for itrace = 1:length(CommonEventGather)
        d_z(:,itrace) = taper(d_z(:,itrace), 5, 5);
        d_z(:,itrace) = bandpassSeis(d_z(:,itrace), dt, lows, highs, 3);

        d_r(:,itrace) = taper(d_r(:,itrace), 5, 5);
        d_r(:,itrace) = bandpassSeis(d_r(:,itrace), dt, lows, highs, 3);
    end

    itmax   = 400;
    minderr = 0.01;
    gauss = 2.5;
    ph      = 5;   % phase delay
    VB      = 0;   % verbose output
    [nt,nx]      = size(d_z);

    itr = zeros(nt, nx);
    for itrace = 1:nx
        R  = d_r(:, itrace);
        Z  = d_z(:, itrace);
        [itr(:, itrace), ~] = makeRFitdecon_la_norm(R, Z, dt, nt, ph, gauss, itmax, minderr);
        CommonEventGather(itrace).RF.itr = itr(:, itrace);
    end

    dist = [];
    for n = 1:length(CommonEventGather)
        dist(end+1) = CommonEventGather(n).TravelInfo.distDeg*(2*pi*6371/360);
    end

    dp = 0.003;
    h  = dist-min(dist);
    nh = length(h);
    p  = -0.05:dp:0.05;
    np = length(p);

    Param.h  = h;
    Param.v  = 1./p;
    Param.nt = nt;
    Param.dt = dt;
    Param.type = 1;
    ma       = zeros(nt, np);

    % Z component
    [~, nx] = size(d_z);

    N1 = 30;   % CG iterations
    N2 = 1;    % Update weight for sparse solution
    [mi_z, ~]  = yc_pcg(@radon_op, Param, d_z, zeros(size(ma)), N1, N2, 1);
    dp_z       = radon_op(mi_z, Param, 1);

    % Plotting
    cmax = 3 * rms(dp_z(:));
    figure('Position', [100 100 1200 800], 'Color', 'w');

    subplot(2,3,1);
    imagesc(h, t, d_z); colormap(seismic(1)); caxis([-cmax cmax]);
    xlabel('Distance (km)'); ylabel('Time (sec)');
    title('Raw Z'); set(gca, 'FontSize', 12);

    subplot(2,3,2);
    imagesc(h, t, dp_z); colormap(seismic(1)); caxis([-cmax cmax]);
    xlabel('Distance (km)'); ylabel('Time (sec)');
    title('Radon Z'); set(gca, 'FontSize', 12);

    subplot(2,3,3);
    imagesc(h, t, d_z - dp_z); colormap(seismic(1)); caxis([-cmax cmax]);
    xlabel('Distance (km)'); ylabel('Time (sec)');
    title('Removed noise (Z)'); set(gca, 'FontSize', 12);

    % R component
    [mi_r, ~] = yc_pcg(@radon_op, Param, d_r, zeros(size(ma)), N1, N2, 1);
    dp_r      = radon_op(mi_r, Param, 1);

    cmax = 3 * rms(dp_r(:));

    subplot(2,3,4);
    imagesc(h, t, d_r); colormap(seismic(1)); caxis([-cmax cmax]);
    xlabel('Distance (km)'); ylabel('Time (sec)');
    title('Raw R'); set(gca, 'FontSize', 12);

    subplot(2,3,5);
    imagesc(h, t, dp_r); colormap(seismic(1)); caxis([-cmax cmax]);
    xlabel('Distance (km)'); ylabel('Time (sec)');
    title('Radon R'); set(gca, 'FontSize', 12);

    subplot(2,3,6);
    imagesc(h, t, d_r - dp_r); colormap(seismic(1)); caxis([-cmax cmax]);
    xlabel('Distance (km)'); ylabel('Time (sec)');
    title('Removed noise (R)'); set(gca, 'FontSize', 12);

    itr_rt = zeros(nt, nx);
    for itrace = 1:nx
        R  = dp_r(1:nt, itrace);
        Z  = dp_z(1:nt, itrace);
        [itr_rt(:, itrace), ~] = makeRFitdecon_la_norm(R, Z, dt, nt, ph, gauss, itmax, minderr);
        CommonEventGather(itrace).RF.itr_rt = itr_rt(:, itrace);
    end

    figure;
    set(gcf,'Position',[0 0 1500 800],'color','w')
    subplot(231)
    wigb(d_z,1,h,t)
    xlabel('Distance')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(a)','Units','normalized','FontSize',18)

    subplot(232)
    wigb(d_r,1,h,t)
    xlabel('Distance')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(b)','Units','normalized','FontSize',18)

    subplot(233)
    %     wigb([sub.itr],4,h,t)
    imagesc(h,t,itr)
    caxis([-0.1 0.1])
    colormap(seismic(1))
    ylim([-2 15])
    xlabel('Distance')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(c)','Units','normalized','FontSize',18)

    subplot(234)
    wigb(dp_z,1,h,t)
    xlabel('Distance')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(d)','Units','normalized','FontSize',18)

    subplot(235)
    wigb(dp_r,1,h,t)
    xlabel('Distance')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(e)','Units','normalized','FontSize',18)

    subplot(236)
    %     wigb([sub.itr1],4,h,t)
    imagesc(h,t,itr_rt)
    caxis([-0.1 0.1])
    colormap(seismic(1))
    ylim([-2 15])
    xlabel('Distance')
    ylabel('Time (sec)')
    set(gca,'fontsize',14)
    text(-0.15,0.98,'(f)','Units','normalized','FontSize',18)

    close all;
end
% Gather results for current event

end
