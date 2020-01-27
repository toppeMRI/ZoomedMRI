function [durations] = setup_duration_sinc(freqs,xfov,yfov,pulse_params)

dt = pulse_params.dt;         % 4usec
gradient_info.slew_rate = pulse_params.slew_rate;      % g/cm/sec, max slew rate
gradient_info.g_max = pulse_params.g_max;              % 4g, maximum gradient amp.
k_max = pulse_params.kz_area;
gam=26751;                   %rad/sec/g
gambar = gam/2/pi;           % Hz/g
up_factor = 100;
dt_upsample = dt/up_factor;
durations = zeros(size(freqs.kx));

g1 = gen_trapezoid(k_max,gradient_info,dt);
g2 = gen_trapezoid(-k_max/2,gradient_info,dt);
a_kx = freqs.kx/xfov;
a_ky = freqs.ky/yfov;

kx_jump_max = 15/xfov; % Hao: original value is 18/xfov
ky_jump_max = 15/yfov;

gx_area = kx_jump_max;
gy_area = ky_jump_max;
gx_cur = gen_trapezoid(gx_area,gradient_info,dt);
gy_cur = gen_trapezoid(gy_area,gradient_info,dt);
pe_time = max(length(gx_cur),length(gy_cur))*dt;

g1_time = length(g1)*dt;
g2_time = length(g2)*dt;
durations(end)=g2_time+pe_time+g1_time/2;
%durations(end)=g2_time+pe_time;

for i=length(durations)-1:-1:1
    %gx_area = (-a_kx(i)+a_kx(i+1));
    %gy_area = (-a_ky(i)+a_ky(i+1));
    durations(i)=durations(i+1)+g1_time+pe_time;
end

return;