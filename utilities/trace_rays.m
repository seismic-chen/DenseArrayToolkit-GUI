function [xray,tray] = trace_rays(raycode,m,p,varargin)
% function to trace the direct wave and later converted wave and calcuate 
% the travel time difference relative to P wave. Code is modified after CREWS
% open source ray tracing code.
% Yunfeng Chen, Dec. 3rd, 2015, Global Seismology Group, University of
% Alberta.
% Input: 
% raycode: two-coloum matrix defines the conversion depth (1st column)
%          and the types of the ray segment (2nd column).
% m      : four-coloum velocity model that have structrue as  
%          depth (1st), vp (2nd), Vs (3rd), Density (4th), the last column
%          could be ignored
% p      : 1*N array defines ray paramert for N seismic wave
% Output:
% Feb. 24, 2017, Yunfeng Chen, add parameter to plot the ray path
% Mar. 22, 2018, Y.C., improve the code to plot PsPms phase
if nargin == 4
    plot_flag = varargin{1};
else
    plot_flag = 0;
end
if nargin == 5
    plot_flag = varargin{1};
    nray = varargin{2};
else
    nray = 1;
end
R = 6371;
zray = raycode(:,1);
raytype=raycode(:,2);
zsrc=zray(1);
zr=zray(end);
zd=max(zray);
%create one-way model
vmod=[];
zmod=[];
zmax=0;

zp = m(:,1);
zmoho = m(end,1);
% scale vs to get vp velocity model
vp = m(:,2);
vs = m(:,3);
% if use RF model or crust 1.0 model, scale vp to get vs
% instead
% vs = m(:,3);
% vp = vs.*kappa(ii);
zs = zp;
vp=vp(:);vs=vs(:);
zp=zp(:);zs=zs(:);
% flat earth transformation
% vpf = (vp*R)./(R-zp);
% zpf = R.*log10(R./(R-zp))./log10(exp(1));
% vsf = (vs*R)./(R-zs);
% zsf = zpf;
% vp=vpf(:);vs=vsf(:);
% zp=zpf(:);zs=zsf(:);

vp=[vp;vp(end)];
zp=[zp;1000000];
zpthk=diff(zp);%thicknesses
vs=[vs;vs(end)];
zs=[zs;1000000];
zsthk=diff(zs);%thicknesses

nseg = 0;
pseg = [];
sseg = [];
for k=1:length(zray)-1
    zr1=zray(k);
    zr2=zray(k+1);
    dir=sign(zr2-zr1);
    incd=round((1+dir)/2);
    incu=round((1-dir)/2);
    if(raytype(k)==1)
        ind1=surround(zp,zray(k));
        ind2=surround(zp,zray(k+1));
        vmod=[vmod;vp(ind1);vp(ind1+dir:dir:ind2)];
        zmod=[zmod;zr1;zp(ind1+incd:dir:ind2+incu)];
        nseg = nseg + length([zr1;zp(ind1+incd:dir:ind2+incu)]);
        pseg = [pseg (nseg - length([zr1;zp(ind1+incd:dir:ind2+incu)]) + 1):nseg];
    else
        ind1=surround(zs,zray(k));
        ind2=surround(zs,zray(k+1));
        vmod=[vmod;vs(ind1);vs(ind1+dir:dir:ind2)];
        zmod=[zmod;zr1;zs(ind1+incd:dir:ind2+incu)];
        nseg = nseg + length([zr1;zs(ind1+incd:dir:ind2+incu)]);
        sseg = [sseg (nseg - length([zr1;zs(ind1+incd:dir:ind2+incu)]) + 1):nseg];
    end
end
zmod=[zmod;zr2];
% zmod=[0;cumsum(abs(diff(zmod)))];
%% ray tracing
xray = [];
tray = [];
[xray,tray] = shootray_corr_mantle(vmod,zmod,p);
% check the bottoming depth of the ray path
if plot_flag == 1
    tt(1,:) = tray(end,:);
    xx(1,:) = xray(end,:);
    %% plot the ray path
    subplot(122)
    xray = [zeros(1,size(xray,2));xray];
    xray_shift = xray - ones(size(xray,1),1)*xx(1,:);
    for l = 1:nseg
        if any(l == pseg)
            if size(xray_shift,2)>=10 % in case there is too many rays, only plot the first 10
                Hp=plot(xray_shift(l:l+1,1:10),zmod(l:l+1),'-k');
            else
                Hp=plot(xray_shift(l:l+1,:),zmod(l:l+1),'-k');
            end
            hold on;
        else
            if size(xray_shift,2)>=10 % in case there is too many rays, only plot the first 10
                Hs=plot(xray_shift(l:l+1,1:10),zmod(l:l+1),'-r');
            else
                Hs=plot(xray_shift(l:l+1,:),zmod(l:l+1),'-r');
            end
            hold on;
        end
    end
    if nray == 5
        set(gcf,'Position',[500 500 1200 500]);
        xlim([-100 5]);
        ylim([0 zmoho]);
        xlabel('Distance (km)');
        ylabel('Depth (km)');
        legend([Hs(1) Hp(1)],'S wave','P wave','Location','NorthWest')
        axis ij;
        set(gca,'FontSize',18);
        subplot(121)
        plot_velocity_model(m(:,1),m(:,2),m(:,3),'thick');
        ylim([0 zmoho]);
        xlim([2 12]);
        legend('S velocity','P velocity');
    end
    
end


function [x,t]=shootray_corr_mantle(v,z,p)
%  Modified after CREW's shootray code
% Mar. 18, 2016, Yunfeng Chen
% Jun. 6, 2016: Note the equation cumsum(thk.*cs./vprop) has alreadly
% considered the mantle correction, so no need to include mantle correction
% term in the calculcation seperately.

%check for critical angle
iprop=1:length(z)-1;
sn = v(iprop)*p;
 %
 % sn is an n by m matrix where n is the length of iprop (the number of
 %    layers propagated through) and m is the length of p (the number of
 %    unique ray parameters to use). Each column of sn corresponds to a
 %    single ray parameter and contains the sin of the vertical angle;
 %
 
 [ichk,pchk]=find(sn>1);

 %compute x and t
 cs=sqrt(1-sn.*sn)+eps;
 vprop=v(iprop)*ones(1,length(p));
 thk=abs(diff(z))*ones(1,length(p));
 if(size(sn,1)>1)
  	x=cumsum( (thk.*sn)./cs);
%     t=cumsum(thk./(vprop.*cs)); 
  	t=cumsum(thk.*cs./vprop); 
%     x=sum( (thk.*sn)./cs);
%  	t=sum(thk./(vprop.*cs));
 else
 	x=(thk.*sn)./cs;
 	t=thk./(vprop.*cs);
 end
 %assign infs
 if(~isempty(ichk))
 	x(pchk)=inf*ones(size(pchk));
	t(pchk)=inf*ones(size(pchk));
 end
