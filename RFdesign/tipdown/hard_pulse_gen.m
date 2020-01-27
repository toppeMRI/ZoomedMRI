function [b1] = hard_pulse_gen(tip_angle,dt)

gam=26751;                          % rad/sec/g
% gambar = gam/2/pi;                  % Hz/g
Nrf = 15; % sample points in rf pulse.
trf = Nrf*dt; 
b1amp = tip_angle/180*pi/gam/trf;  %gauss
b1 = ones(Nrf,1)*b1amp; 
disp(sprintf('unit pulse length %f ms',length(b1)*dt*10^3));
return;