function leastSquaresMig(DataStruct,param)
addpath utils/
addpath ssa/
addpath radon/
addpath OneWay_RF/
addpath export_fig/
%% load paramerters
is_radon=param.is_radon;
is_ssa=param.is_ssa;

% location of profile
lon1 = param.lon1;
lat1 = param.lat1;
lon2 = param.lon2;
lat2 = param.lat2;

% grid info
dx = param.dx;
nx = param.nx;
xpad=param.xpad;
xmax=param.xmax;

minsnr = param.minsnr;
mintrace = param.mintrace;
%% prepare velocity model for migration
% [vel,vel_s,x,z] = rflsm_create_initial_model(param);

% find qualified events
slat = cellfun(@(stationinfo) stationinfo.stla, {DataStruct.StationInfo}, 'UniformOutput', false);
slon = cellfun(@(stationinfo) stationinfo.stlo, {DataStruct.StationInfo}, 'UniformOutput', false);
slat = cell2mat(slat)';
slon = cell2mat(slon)';
[~,idx,~]=unique([slat slon],'rows');
slat = slat(idx);
slon = slon(idx);

elat = cellfun(@(eventinfo) eventinfo.evla, {DataStruct.EventInfo}, 'UniformOutput', false);
elon = cellfun(@(eventinfo) eventinfo.evlo, {DataStruct.EventInfo}, 'UniformOutput', false);
eventid = cellfun(@(eventinfo) eventinfo.evid, {DataStruct.EventInfo}, 'UniformOutput', false);
elat = cell2mat(elat)';
elon = cell2mat(elon)';
[~,idx,~]=unique([elat elon],'rows');
elat = elat(idx);
elon = elon(idx);
eventid = eventid(idx);
max_angle_diff = 15;
consistent_earthquakes_index = filter_earthquakes_by_azimuth(slon, slat, elon, elat, max_angle_diff);
eventid = eventid(consistent_earthquakes_index);

% project the stations onto a profile
[line_points, projected_points] = fit_and_project_stations(slon, slat,1);
slatp=projected_points(:,2);
slonp=projected_points(:,1);
lat1 = line_points(1,2);
lon1 = line_points(1,1);
lat2 = line_points(end,2);
lon2 = line_points(end,1);

[deg0,~]= distance(lat1,lon1,slatp,slonp);
rx = deg0*2*pi*6371/360;
%% create velocity model
xpad = 50;
model = obtain_crust1_QB();
dx = 4;
x = 0-xpad:dx:350+xpad;
nx = length(x);
dz = 1;
z = 0:dz:100;
nz = length(z);
vp = interp1(model(:,1),model(:,2),z,'nearest','extrap');
vs = interp1(model(:,1),model(:,3),z,'nearest','extrap');
vel = repmat(vp(:),1,nx);
vel_s = repmat(vs(:),1,nx);

% smooth the velocity model
N=5;
[vel,~]=moving_avg(vel,N,'constant',2);
[vel,~]=moving_avg(vel,N,'constant');
[vel_s,~]=moving_avg(vel_s,N,'constant',2);
[vel_s,~]=moving_avg(vel_s,N,'constant');

figure;
imagesc(x,z,vel); hold on; scatter(rx,zeros(size(rx)),100,'r^','filled')
%% loop over all shots (events)
% eventid=getEvents(DataStruct);
nshot=length(eventid);
dmig = zeros(nz,nx,nshot);
dmigls = zeros(nz,nx,nshot);
empty_shot_index = [];
preserved_shot_id = [];
for ishot=1:nshot
    disp(['Processing shot ', num2str(ishot)]);
    gather = getCommonEventGather(DataStruct,eventid{ishot});
    if length(gather)<50
        empty_shot_index = [empty_shot_index ishot];
        continue
    end
    preserved_shot_id = [preserved_shot_id; eventid{ishot}];
    itr = cellfun(@(rf) rf.itr, {gather.RF}, 'UniformOutput', false);
    itr = cell2mat(itr);

    slat = cellfun(@(stationinfo) stationinfo.stla, {gather.StationInfo}, 'UniformOutput', false);
    slon = cellfun(@(stationinfo) stationinfo.stlo, {gather.StationInfo}, 'UniformOutput', false);
    slat = cell2mat(slat)';
    slon = cell2mat(slon)';

    % project the stations onto a profile
    [line_points, projected_points] = fit_and_project_stations(slon, slat);
    slatp=projected_points(:,2);
    slonp=projected_points(:,1);
    % calculate the receiver location along the profile
    [deg0,~]= distance(lat1,lon1,slatp,slonp);
    rx = deg0*2*pi*6371/360;

    param.rx = rx;

    % remove bad traces
%     remove= [gather.RF] < minsnr;
%     gather(remove)=[];
%     if length(gather) < mintrace
%         continue
%     end
    %% Binning
    nt = length(gather(1).TimeAxis.t_resample);
    d = zeros(nt,nx);
    for i=1:nx
        idx=rx>=x(i)-dx/2 & rx<=x(i)+dx/2;
        if sum(idx)>0
            d(:,i) = mean(itr(:,idx),2);
        end
    end
    figure;
    wigb(d,2,x,gather(1).RF.ittime)
    ylim([-2 30])
    %% Parameters for SSA/CAZDOW reconstruction
    param.flow = 0.1;    % minimum frequency
    param.fhigh = 1.2;  % maximum frequency
    param.rank_p = 8;  % rank value

    param.alpha=0.9;    % trade-off parameter
    param.n_iter=20;    % iteration number

    param.xmax = max(x);

    ittime = gather(1).RF.ittime;
    dt = gather(1).TimeAxis.dt_resample;
    
    idx = ittime < 60;
    d = d(idx,:);
    ittime = ittime(idx);

    is_ssa = 0;
    if is_ssa
        i1=xpad/dx+1;
        i2=nx-xpad/dx-1;
        din = d(:,i1:i2);
        [dout] = rflsm_ssa(din,dt,param);
        d(:,i1:i2) = dout;
    end
    %% forward propagation of P wavefield
    vp = mean(vel(end,:));
    raypAll = cellfun(@(travelinfo) travelinfo.rayParam/6371, {gather.TravelInfo}, 'UniformOutput', false);
    bazAll =  cellfun(@(travelinfo) travelinfo.baz, {gather.TravelInfo}, 'UniformOutput', false);
    rayp = mean(cell2mat(raypAll));
    baz = mean(cell2mat(bazAll));

    [nt,nx] = size(d);
    param.gauss = 2.5;
    param.ph = 5;
    param.nx = nx;
    param.dx = dx;
    
    [src,pos,tshift] = rflsm_create_src(dt,nt,rayp,baz,vp,param);
 
    % taper RF to remove later conversions
%     wl=floor(((40-ittime(1))/dt+1));
%     w = [tukeywin(wl,0.75); zeros(size(d,1)-wl,1)];
%     win=w*ones(1,size(d,2));
%     d=d.*win;

    % scale the src
%     src = src *max(d(:));
    src = src*max(mean(itr,2));

    param.dx = dx;
    param.dz = dz;
    param.x = x;
    param.z = z;
    param.nx = nx;
    param.nz = nz;
    param.xmax = max(x);
    param.flow = 0.1;
    param.fhigh = 1.2;

    dshift = rflsm_shift_rfs(d,ittime,vel,vel_s,src,pos,tshift,param);
    %% Migration
    [mig,~] = rflsm_migration(dshift,ittime,vel,vel_s,src,pos,tshift,param);
    dmig(:,:,ishot) = mig;
    %% LSM
    param.itermax = 20;
    param.mu = 0.1;
    [migls,~] = rflsm_lsm(dshift,ittime,vel,vel_s,src,pos,tshift,param);
    dmigls(:,:,ishot) = migls;
    %% close all figures
     close all
end
dmig(:,:,empty_shot_index) = [];
dmigls(:,:,empty_shot_index) = [];
%% save migration results
% save(fullfile(figdir,'mig.mat'), 'dmig', 'x', 'z')
% save(fullfile(figdir,'migcg.mat'), 'dmigls', 'x', 'z')
d2d=mean(dmig,3);
d2dls=mean(dmigls,3);

%     ngrid_x = 4;
%     ngrid_y = 4;
%     kernel = ones(ngrid_x,ngrid_y)/(ngrid_x * ngrid_y);
%     d2d = conv2(d2d, kernel, 'same');
%     d2dls = conv2(d2dls, kernel, 'same');
%% plot migration results
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
cmax=0.2;
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
colormap(seismic(3));
text(-0.12,0.98,'b)','Units','normalized','FontSize',18)
figname='mig_compare_real.png';
% export_fig(fullfile(figdir,figname));