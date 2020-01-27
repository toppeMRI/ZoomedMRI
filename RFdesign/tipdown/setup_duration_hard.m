function [durations] = setup_duration_hard(freqs,xfov,yfov,zfov,pulse_params)

dt = pulse_params.dt;         % 4usec
gradient_info.slew_rate = pulse_params.slew_rate;      % g/cm/sec, max slew rate
gradient_info.g_max = pulse_params.g_max;              % 4g, maximum gradient amp.
% k_max = pulse_params.kz_area;
gam=26751;                   %rad/sec/g
gambar = gam/2/pi;           % Hz/g


% duration of phase encoding part: 
kx_jump_max = pulse_params.max_pe_jump; % Hao: original value is 18/xfov
ky_jump_max = kx_jump_max;
kz_jump_max = kx_jump_max;
gx_area = kx_jump_max;
gy_area = ky_jump_max;
gz_area = kz_jump_max;
gx_cur = gen_trapezoid(gx_area,gradient_info,dt);
gy_cur = gen_trapezoid(gy_area,gradient_info,dt);
gz_cur = gen_trapezoid(gz_area,gradient_info,dt);
pe_time = max(max(length(gx_cur),length(gy_cur)),length(gz_cur))*dt;

% duration of RF pulse part: 
rf_time = 20*dt; % 60us hard pulse

% total duration of each subpulse
durations = zeros(size(freqs.kx));
durations(end)= pe_time + rf_time;
for i=length(durations)-1:-1:1
    durations(i)=durations(i+1)+rf_time+pe_time;
end

return;