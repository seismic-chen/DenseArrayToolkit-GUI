function [HKresults, plotData] = HKstacking(DataStruct, app, Hall, kappa)
    % HKstacking - 利用接收函数Hk叠加计算台站地壳厚度和Vp/Vs(泊松比)
    %
    % 输入:
    %   DataStruct - 包含台站和事件数据的结构体
    %
    % 输出:
    %   HKresults：['Sta','Nrfs','Lat','Lon','Ele','H','Hstd','k','kstd']
    %   H    - 台站下方地壳厚度
    %   Hstd     - 台站下方地壳厚度误差估计
    %   k - 台站下方地壳平均Vp/Vs
    %   kstd - 台站下方地壳平均Vp/Vs误差估计
    global Hall
    global kappa

    R = 6371; % radius of the Earth
    edep = 0;
%     Hall = 30:0.1:70;
    % for figure plot
%     kappa = 1.6:0.01:2.0;
    % structure for HK results
    hkstat = struct();
    if_print_result = 1;
    if if_print_result == 1
        fid1 = fopen('./results/HK_bootstrap_tmp.txt','a');
        fprintf(fid1,'%-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s\n',...
            'Sta','Nrfs','Lat','Lon','Ele','H','Hstd','k','kstd','Hb_mean','Hb_std','kb_mean','kb_std');
    end

    station = getStations(DataStruct); 
    stationList = {station.sta};
    nsta = length(stationList);
    for n = 1:nsta
        gather = getCommonStationGather(DataStruct, stationList{n});
        if isempty(gather)
            continue;  % 跳过空数据
        end
    end

    % 初始化结构体数组
    HKresults = struct('Sta', [],'Nrfs', [],'Lat', [],'Lon', [],'Ele', [],'H', [],'Hstd', [],'k', [],'kstd', []);
%     plotData = struct('Hall', {}, 'kappa', {}, 'Cn', {}, 'besth', {}, 'bestk', {}, 'hstd', {}, 'kstd', {});  % 用于保存绘图数据的结构体
    plotData = repmat(struct('Hall',[],'kappa',[],'Cn',[],'besth',[],'bestk',[],'hstd',[],'kstd',[],...
                              't',[],'seis',[],'dist',[],'rayp',[],'Tdelay_Ps',[], ...
                              'Tdelay_PpPms',[],'Tdelay_PpSms',[], 'vp_ave',[], 'indexh',[], 'indexk',[]), 1, n);%indexh,indexk
    for n = 1:length(stationList)
        % 允许 UI 更新，检测用户是否按下停止按钮
        drawnow;  
        if app.StopFlag
            % 如果停止标志为 true，则退出循环
            disp(['计算在第 ', num2str(n), ' 个台站时被停止。']);
            break;
        end
        gather = getCommonStationGather(DataStruct, stationList{n});
        if isempty(gather)
            continue;  % 跳过空数据
        end
        seis=[];rayp=[];dist=[];
        net = gather(1).StationInfo.network;
        sta = gather(1).StationInfo.sta;
        lat = gather(1).StationInfo.stla;
        lon = gather(1).StationInfo.stlo;
        ele = gather(1).StationInfo.stel;
        hkstat.station(n).slat = lat;
        hkstat.station(n).slon = lon;
        hkstat.station(n).sele = ele;
        hkstat.station(n).num_rf = length(DataStruct(n).RF);

        for nrf = 1:length(gather)
            seis = [seis gather(nrf).RF.itr];
        end
        for nrf = 1:length(gather)
            rayp = [rayp gather(nrf).TravelInfo.rayParam];
        end
        for nrf = 1:length(gather)
            dist = [dist gather(nrf).TravelInfo.distDeg];
        end

        t = gather(1).RF.ittime;
        dt = t(2)-t(1);
%         dist = gather(n).TravelInfo.distDeg;
%         rayp = gather(n).TravelInfo.rayParam;
        Pflat = rayp/6371;
%         Pflat = rayp_to_dist(rayp','km');

%         seis = gather(n).RF.itr;
        keep = any(seis);
        seis = seis(:,keep);
        rayp = rayp(keep);
        dist = dist(keep);
        % initialize some variable
        Tdelay_Ps = zeros(length(Hall),length(Pflat),length(kappa));
        Tdelay_PpPms = zeros(size(Tdelay_Ps));
        Tdelay_PpSms = zeros(size(Tdelay_Ps));
        stack1 = zeros(length(Hall),length(Pflat));
        stack2 = zeros(length(Hall),length(Pflat));
        stack3 = zeros(length(Hall),length(Pflat));
        %% search for H
        m0 = [];
        nsedi = 0;
        [m0,nsedi] = obtain_crust1(lat,lon);
        % nsedi indicates the number of sedimentary layer,e.g., nsedi = 2 means
        % the top two layers are sediments
        dmoho = m0(end,1);
        vp_mantle = m0(end,2);
        vs_mantle = m0(end,3);
        rho_mantle = m0(end,4);
        tic;
        for i = 1:length(Hall)
            m = [];
            disp(['Now search H = ',num2str(Hall(i))]);
            H = Hall(i);
            % find the index of depth that is closest to H
            index_H = surround(m0(:,1),H);
            if isempty(index_H)
                index_H = size(m0,1);
            end
            m = m0(1:index_H,:); % only keep the model above H
            % if H == dmoho, then do not need to add this layer
            if H ~= dmoho
                m(end+1,1:4) = [H vp_mantle vs_mantle rho_mantle];
            end
            if H > dmoho % set the velcoity of layer above H to lower crust velocity
                m(end-1,2:4) = m(end-2,2:4);
            end
    
            %% define the velocity model
    %         m = obtain_jeff_model(H,lat,lon);
            %% define the ray code to trace the converted wave and multiples
            raycode_P = [H-eps*10000 1; 0+eps*10000 2];
            raycode_Pms = [H-eps*10000 2; 0+eps*10000 2];
            raycode_PpPms = [H-eps*10000 1; 0+eps*10000 1; H-eps*10000 2; 0+eps*10000 2];
            raycode_PpSms = [H-eps*10000 1; 0+eps*10000 2; H-eps*10000 2; 0+eps*10000 2];
            raycode_all{1} = raycode_P;
            raycode_all{2} = raycode_Pms;
            raycode_all{3} = raycode_PpPms;
            raycode_all{4} = raycode_PpSms;
            tt = zeros(4,length(Pflat));
            xx = zeros(4,length(Pflat));
            for ii = 1:length(kappa)
                mtmp = m;
                % do not change the velocities in the mantle and sediments
                mtmp(nsedi+1:end-1,3) = mtmp(nsedi+1:end-1,2)/kappa(ii);
                for nray = 1:length(raycode_all)
                    raycode = raycode_all{nray};
                    % trace the rays
                    [xtmp,ttmp] = trace_rays(raycode,mtmp,Pflat);
                    tt(nray,:) = ttmp(end,:);
                    xx(nray,:) = xtmp(end,:);
                end
                %% calucalte the corrected time
                % tt is 4 by N array contains the travel time of 1 P and 3 S phases 
                % from N RFs
	            Tdelay = zeros(1,3);
                for j = 1:size(seis,2)
                    Tdelay = [tt(2,j)-tt(1,j); tt(3,j)-tt(1,j); tt(4,j)-tt(1,j)];
                    Tdelay_Ps(i,j,ii) = Tdelay(1);
                    Tdelay_PpPms(i,j,ii) = Tdelay(2);
                    Tdelay_PpSms(i,j,ii) = Tdelay(3);
                    % find the amplitude in the seismogram
                    x = seis(:,j);
                    % Ps timing
                    tc1 = abs(t(1)) + Tdelay(1);
                    nc1 = round((tc1)/dt);
                    nmin1 = round((tc1-0.1)/dt);
                    nmax1 = round((tc1+0.1)/dt);
                    % take the amplitue
                    stack1(i,j,ii) =  mean(x(nmin1:nmax1));
                    % PpPms timing
                    tc2 = abs(t(1)) + Tdelay(2);
                    nc2 = round((tc2)/dt);
                    nmin2 = round((tc2-0.1)/dt);
                    nmax2 = round((tc2+0.1)/dt);
                    % take the amplitue
                    stack2(i,j,ii) =  mean(x(nmin2:nmax2));
                    % PpSms timing
                    tc3 = abs(t(1)) + Tdelay(3);
                    nc3 = round((tc3)/dt);
                    nmin3 = round((tc3-0.1)/dt);
                    nmax3 = round((tc3+0.1)/dt);
                    % take the amplitue, negative so need '-' sign
                    stack3(i,j,ii) =  -mean(x(nmin3:nmax3));
                end
            end
            toc;
        end
        %% sum up the row vector
        % Nth root stacking, N = 1 is the linear stack
        K = size(seis,2);
        N = 1;
        Cn1 = zeros(length(Hall),length(kappa));
        Cn2 = zeros(length(Hall),length(kappa));
        Cn3 = zeros(length(Hall),length(kappa));
        rd1 = zeros(1,length(Hall));
        rd2 = zeros(1,length(Hall));
        rd3 = zeros(1,length(Hall));
        for ii = 1:size(kappa,2)
            for i = 1:size(Hall,2)
                % Pms
                rd1(i) = 1/K*sum(sign(stack1(i,1:K,ii)).*(abs(stack1(i,1:K,ii))).^(1/N));
                Cn1(i,ii) = rd1(i)*abs(rd1(i))^(N-1);
                % PpPms
                rd2(i) = 1/K*sum(sign(stack2(i,1:K,ii)).*(abs(stack2(i,1:K,ii))).^(1/N));
                Cn2(i,ii) = rd2(i)*abs(rd2(i))^(N-1);
                % PpSms
                rd3(i) = 1/K*sum(sign(stack3(i,1:K,ii)).*(abs(stack3(i,1:K,ii))).^(1/N));
                Cn3(i,ii) = rd3(i)*abs(rd3(i))^(N-1);
            end
        end
        vp_ave = m0(end,1)/sum(diff(m0(:,1))./m0(1:end-1,2));
        %% stacking the results for 3 phases
        %w1 = 0.6;
        %w2 = 0.3;
        %w3 = 0.1;
        w1 = 1/3;
        w2 = 1/3;
        w3 = 1/3;
        Cn = [w1.*Cn1+w2.*Cn2+w3.*Cn3];
        Cn = Cn/max(Cn(:));
        Cn(Cn<0) = 0;
        [~,index] = max(Cn(:));
        [indexh,indexk] = ind2sub(size(Cn),index);
        besth = Hall(indexh);
        bestk = kappa(indexk);
        % save to hkstat structure
        hkstat.station(n).H_raytracing = besth;
        hkstat.station(n).k_raytracing = bestk;    
        
        % calcualte the uncertainty using bootstraping method
        nboot = 200; %number of trial
        % To use the matlab built-in function, we need to convert HK stacking
        % results for each phase, which are NH*N*NK matrix to a 2D array with
        % dimension N by NH*NK such that each row represents the HK results of a
        % trace
        [NH,N,NK] = size(stack1);
        tmp1 = reshape(permute(stack1,[1,3,2]),NH*NK,N);
        tmp2 = reshape(permute(stack2,[1,3,2]),NH*NK,N);
        tmp3 = reshape(permute(stack3,[1,3,2]),NH*NK,N);
        % tranpse the matrix to N by NH*NK
        tmp1 = tmp1';
        tmp2 = tmp2';
        tmp3 = tmp3';
        % apply bootstrap
        [bootstat,bootsam] = bootstrp(nboot,@hk_stack_bootstrap,tmp1,tmp2,tmp3);
        hstd = std(bootstat(:,1));
        kstd = std(bootstat(:,2));  
        % remove the 1/3 of the values that are far away from the meadian picked
        H_bootstrap = bootstat(:,1);
        k_bootstrap = bootstat(:,2);
        H_mean = mean(H_bootstrap);
        k_mean = mean(k_bootstrap);
        H_median = median(H_bootstrap);
        k_median = median(k_bootstrap);
        % save to hkstat structure
        hkstat.station(n).H_std = hstd;
        hkstat.station(n).k_std = kstd;
        hkstat.station(n).H_mean = H_mean;
        hkstat.station(n).k_mean = k_mean;
        hkstat.station(n).H_median = H_median;
        hkstat.station(n).k_median = k_median;
        
        % calculate the difference between each bootstrap value and the median
        H_diff = abs(H_bootstrap-H_median);
        [~,order] = sort(H_diff,'descend');
        index_start = round(length(order)/3,0);
        index_keep = order(index_start:end);
        % index_keep = order(1:end);
        % bootstrap mean
        H_bootstrap_mean = mean(H_bootstrap(index_keep));
        H_bootstrap_median = median(H_bootstrap(index_keep));
        H_bootstrap_std = std(H_bootstrap(index_keep));
        k_bootstrap_mean = mean(k_bootstrap(index_keep));
        k_bootstrap_median = median(k_bootstrap(index_keep));
        k_bootstrap_std = std(k_bootstrap(index_keep));
        hkstat.station(n).H_std_no_outlier = H_bootstrap_std;
        hkstat.station(n).k_std_no_outlier = k_bootstrap_std;
        hkstat.station(n).H_mean_no_outlier = H_bootstrap_mean;
        hkstat.station(n).k_mean_no_outlier = k_bootstrap_mean;
        hkstat.station(n).H_median_no_outlier = H_bootstrap_median;
        hkstat.station(n).k_median_no_outlier = k_bootstrap_median;
    %     H_bootstrap_mean = mean(H_bootstrap);
    %     H_bootstrap_std = std(H_bootstrap);
    %     k_bootstrap_mean = mean(k_bootstrap);
    %     k_bootstrap_std = std(k_bootstrap);
        % print results 
        nrfs = size(seis,2);
        if nrfs < 3
            continue
        end
        nrfs = size(seis,2);
        if if_print_result == 1
            fprintf(fid1,'%-7s %-7d %-7.2f %-7.2f %-7.2f %-7.2f %-7.2f %-7.2f %-7.2f %-7.2f %-7.2f %-7.2f %-7.2f\n',...
                sta,nrfs,lat,lon,ele,besth,hstd,bestk,kstd,H_bootstrap_mean,H_bootstrap_std,k_bootstrap_mean,k_bootstrap_std);
        end
        % 每次计算完当前台站的结果后，将其保存到结构体数组中
        % 'Sta','Nrfs','Lat','Lon','Ele','H','Hstd','k','kstd'
        HKresults(n).Sta = sta;
        HKresults(n).Nrfs = nrfs;
        HKresults(n).Lat = lat;
        HKresults(n).Lon = lon;
        HKresults(n).Ele = ele;
        HKresults(n).H = besth;
        HKresults(n).Hstd = hstd;
        HKresults(n).k = bestk;
        HKresults(n).kstd = kstd;
        %用来在UI中画图

        plotData(n).Hall = Hall;
        plotData(n).kappa = kappa;
        plotData(n).Cn = Cn;
        plotData(n).besth = besth;
        plotData(n).bestk = bestk;
        plotData(n).hstd = hstd;
        plotData(n).kstd = kstd;
        plotData(n).t = t;
        plotData(n).seis = seis;
        plotData(n).dist = dist;
        plotData(n).rayp = rayp;
        plotData(n).Tdelay_Ps = Tdelay_Ps;
        plotData(n).Tdelay_PpPms = Tdelay_PpPms;
        plotData(n).Tdelay_PpSms = Tdelay_PpSms;
        plotData(n).vp_ave = vp_ave;
        plotData(n).indexh = indexh;
        plotData(n).indexk = indexk;
        
%         disp(plotData)
        %% plot
        fig1 = figure('Visible','off');
        set(fig1,'Position',[500 500 1000 1000]);
        imagesc(Hall,kappa,Cn');hold on;

        textinfoh = ['Depth=',num2str(besth,'%4.1f'), setstr(177),num2str(hstd,'%2.1f'),' km'];
        textinfok = [ 'Vp/Vs=',num2str(bestk,'%4.2f'),setstr(177),num2str(kstd,'%4.2f')];
        scatter(besth,bestk,100,'white','x');hold off;
        th = annotation('textbox','String',textinfoh,'Color','k','FontSize',24,'EdgeColor','none','BackgroundColor','w');
        tk = annotation('textbox','String',textinfok,'Color','k','FontSize',24,'EdgeColor','none','BackgroundColor','w');
        
        set(th,'Position',[0.655 0.875 0.25 0.05])
        set(tk,'Position',[0.655 0.825 0.25 0.05])

        axis xy;
        xlabel('Depth(km)','FontSize',24);
        ylabel('Vp/Vs','FontSize',24);
        title(sta,'FontSize',24);
        set(gca,'FontSize',20);
        % shape 0 value to while color
        cmap = colormap(jet);
        %cmap = [1 1 1; cmap];
        colormap(cmap)
        figname1 = ['./figures/HK/',sta,'_HK_raytracing.fig'];
        figname = ['./figures/HK/',sta,'_HK_raytracing'];
        saveas(fig1,figname,'png');
        saveas(fig1,figname,'fig');
        % export_fig(fig2,figname,'-png');

        fig2 = figure('Visible','off');
        set(fig2,'Position',[500 500 1000 1000]);
        for id = length(dist):-1:1
            % plot traces
            x = seis(:,id)*30;
            shift = dist(id);
            p1 = plot(t,x+shift);hold on;
            % fill in the color for positive phase
            thre = 0.01;
            upper = x; upper(upper<=thre) = 0;
            upper(1) = 0; upper(end) = 0;
            lower = zeros(length(x),1);
            jbfill(t(:),upper+shift,lower+shift,'r','k',1,1.0);hold on;
            
            upper = zeros(length(x),1);
            lower = x; lower(lower>=-thre)=0;
            lower(1) = 0; lower(end) = 0;
            jbfill(t(:),upper+shift,lower+shift,[0.17,0.17,0.17],'k',1,1.0);hold on;
            xlim([-5,35]);
            ylim([20,100]);
            set(p1,'LineWidth',1)
        end
        H = besth; k = bestk;
        for id = length(dist):-1:1
            % mark the arrival time of major phases based on the calculated HK
            % results
            shift = dist(id);
            p = rayp(id);
            t0p1s = H * (sqrt( (k/vp_ave)^2 - p^2 ) - sqrt( (1/vp_ave)^2 - p^2 ));
            t2p1s = H * (sqrt( (k/vp_ave)^2 - p^2 ) + sqrt( (1/vp_ave)^2 - p^2 ));
            t1p2s = 2 * H * sqrt( (k/vp_ave)^2 - p^2 );
            scatter(0,shift,40,'blue','filled');
            scatter(t0p1s,shift,40,'black','filled');
            scatter(t2p1s,shift,40,'black','filled');
            scatter(t1p2s,shift,40,'black','filled');
            t0p1s = Tdelay_Ps(indexh,id,indexk);
            t2p1s = Tdelay_PpPms(indexh,id,indexk);
            t1p2s = Tdelay_PpSms(indexh,id,indexk);
            scatter(0,shift,40,'blue','filled');
            scatter(t0p1s,shift,40,'blue','filled');
            scatter(t2p1s,shift,40,'blue','filled');
            scatter(t1p2s,shift,40,'blue','filled');
        end
        hold off;
        title(sta,'FontSize',24);
        xlabel('Time (sec)');
        ylabel('Distance (deg)');
        set(gca,'FontSize',20);
        figname = ['./figures/HK/',sta,'_waveform'];
        figname2 = ['./figures/HK/',sta,'_waveform.fig'];
        saveas(fig2,figname,'fig');
        saveas(fig2,figname,'png');

        close all;
    end
    
    if if_print_result == 1
        fclose(fid1);
    end
end
    