function plotrf(b1,gx,gy,gz)
% function plotrf(b1,gx,gy,gz)
% $Id: plotrf_sos3d.m,v 1.1 2015/05/22 21:39:24 jfnielse Exp $

ncoils = size(b1,2);
nt = size(b1,1);
dt = 4e-3;  % ms
T = linspace(dt/2,nt*dt-dt/2,nt);

cols = 'bbbbbbbb';

figure;
for c = 1:ncoils
	subplot(2,ncoils,c);        plot(T,abs(b1(:,c)),cols(c));   ylabel(['coil ' num2str(c) ', abs(b1) (Gauss)']);
	xlabel('time (msec)');
	subplot(2,ncoils,c+ncoils); plot(T,angle(b1(:,c)),cols(c)); ylabel(['coil ' num2str(c) ', angle(b1)']);
	xlabel('time (msec)');
end

figure;
subplot(1,3,1); plot(T,gx); ylabel('gx (G/cm)');
	xlabel('time (msec)');
subplot(1,3,2); plot(T,gy); ylabel('gy (G/cm)');
	xlabel('time (msec)');
subplot(1,3,3); plot(T,gz); ylabel('gz (G/cm)');
	xlabel('time (msec)');

return;
% EOF
