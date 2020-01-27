function [pulses,gx,gy,gz] = full_pulse_gen_sinc(tip_angle,slice_thickness,kx,ky,kxy_max,kz_max,weights,ncoils,slew_rate,g_max,dt)

gradient_info.slew_rate = slew_rate;
gradient_info.g_max = g_max;

gx = [];
gy = [];
gz = [];
pulses = [];

n_seg = length(kx);
ncoils = size(weights,2);

b1 = sinc_pulse_gen(tip_angle,slice_thickness,kz_max,slew_rate,g_max,dt);

b1 = b1*50/49; %amplitude correction

gz1 = gen_trapezoid(kz_max,gradient_info,dt);
gz2 = gen_trapezoid(-kz_max/2,gradient_info,dt);

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
    
    gz_new = gz1;
    gx = [zeros(size(gz_new));gx_new;gx];
    gy = [zeros(size(gz_new));gy_new;gy];    
    gz = [gz_new;zeros(size(gx_new));gz];
    b_new = [b1;zeros(size(gx_new))];
    pulses_new = [];
    for j = 1:ncoils
        pulses_new = [pulses_new b_new*weights(i,j)];
    end
    pulses = [pulses_new;pulses];
        
    gz1 = -gz1;
end

gz = [gz;gz2];
gx = [gx;zeros(size(gz2))];
gy = [gy;zeros(size(gz2))];
pulses_new = zeros(length(gz2),ncoils);
pulses = [pulses;pulses_new];

return;

