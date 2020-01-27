function [pulses,gx,gy,gz] = full_pulse_gen_hard(tip_angle,slice_thickness,kx,ky,kxy_max,kz_max,weights,ncoils,slew_rate,g_max,dt)

gradient_info.slew_rate = slew_rate;
gradient_info.g_max = g_max;

gx = [];
gy = [];
gz = [];
pulses = [];

n_seg = length(kx);
ncoils = size(weights,2);

%b1 = hard_pulse_gen(tip_angle,dt);

gam = 26751;                          % rad/sec/g
% gambar = gam/2/pi;                  % Hz/g
Nrf = 15; % sample points in one hard pulse.
trf = Nrf*dt; 
b1amp = sin(tip_angle/180*pi)/trf/gam/1e4;  %gauss
b1 = ones(Nrf,1)*b1amp; 
disp(sprintf('unit pulse length %f ms',length(b1)*dt*10^3));

%b1 = b1*50/49; %amplitude correction

for i = n_seg:-1:1
    if ( i == n_seg )
        dkx = kx(end);
        dky = ky(end);
    else
        dkx = kx(i) - kx(i+1);
        dky = ky(i) - ky(i+1);        
    end

    
    gx_new = gen_trapezoid(-dkx,gradient_info,dt);
    gy_new = gen_trapezoid(-dky,gradient_info,dt);
    
    %g_xy_len = max(length(gx_new),length(gy_new));
    g_xy_len = length(gen_trapezoid(kxy_max,gradient_info,dt));
    gx_new = [gx_new;zeros(g_xy_len-length(gx_new),1)];
    gy_new = [gy_new;zeros(g_xy_len-length(gy_new),1)];
    

    gx = [zeros(size(b1));gx_new;gx];
    gy = [zeros(size(b1));gy_new;gy];    

    b_new = [b1;zeros(size(gx_new))]; % pulse = weight/gam/Trf*sin(alpha)
    pulses_new = [];
    for j = 1:ncoils
        pulses_new = [pulses_new b_new*weights(i,j)];
    end
    pulses = [pulses_new;pulses];
        

end


% gx = [gx;zeros(size(gz2))];
% gy = [gy;zeros(size(gz2))];
gz = zeros(size(gx)); 
% pulses_new = zeros(length(gz2),ncoils);
% pulses = [pulses;pulses_new];

return;

