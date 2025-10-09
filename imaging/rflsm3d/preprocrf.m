function [itr,src_func] = preprocrf(rf0,param)

TIME = param.Ti;
t1 = -3;
t2 = 20;
dt = param.dt;
nt1 = abs(t1)/dt;
nt2 = t2/dt;
x = size(rf0,2);

% normaliztion
itr = rf0./max(rf0(:));
nfft = size(itr,1);

% extract wavelet
% srctmp=mean(itr,2);
% [win] = waveform_win(srctmp,TIME,t1,1,-1);
% src=srctmp.*win;
% src_func = src./max(src(:));

gauss = param.gauss;
ph = param.phaseshift;
gaussF = gaussFilter( dt, nfft, gauss );
% phase shift in radians
shift_i = round(ph/dt); % removed +1 from here.
p = 2*pi*(1:nfft).*shift_i./(nfft);

% apply shift
Xf = gaussF.*(cos(p) - 1i .* sin(p) );

% back into time
src = real( ifft(Xf, nfft) )/cos(2*pi*shift_i/nfft);

delf = (1/dt)/nfft;
gnorm = sum(gaussF)*delf*dt;
src = real(src(:))/gnorm;
src_func = src./max(src(:));

% taper RF to remove later conversions
[win] = waveform_win(src_func,TIME,t1,t2-5,3);
win = win*ones(1,size(itr,2));
itr = itr.*win;

% src_func = gradient(src_func,-0.5);


end