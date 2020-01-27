function grad_waveform = gen_trapezoid2(k, gradient_info, dt)
%k : cycle/cm
%slew_rate : g/cm/sec
%g_max : g/cm
%dt : sec

slew_rate = gradient_info.slew_rate;
g_max = gradient_info.g_max;

gam=26751;                   %rad/sec/g
gambar = gam/2/pi;           % Hz/g

if ( k == 0 )
    grad_waveform = [];
    return;
end

rise_time = g_max/slew_rate;
ramp = [0:1:ceil(rise_time/dt)]*slew_rate*dt;
ramp(end)=g_max;
threshold_k = sum(ramp)*dt*gambar*2;

sign_k = sign(k);
k = abs(k);

if ( k <= threshold_k )    
    %rise_time = sqrt(k/gambar/slew_rate);
    N1 = ceil(sqrt(k/gambar/(dt^2)/slew_rate));
    for n = N1:-1:1
        est_k = n*(n+1)*(dt^2)*slew_rate*gambar;
        if (est_k < k)
            N1 = n;
            break;
        end
    end
    time = [0:1:N1]*dt;
    rise_waveform = time*slew_rate;
    if (rise_waveform(end) >= g_max)
        rise_waveform(end)=g_max;
    end

    rise_sum = sum(rise_waveform)*dt*gambar;
    k_err = k/2 - rise_sum;

    for l = 1:length(rise_waveform)
        amp_needed = k_err/(l*dt)/gambar;
        idx = find(amp_needed<rise_waveform);
        if ( length(idx) == 0 )
            continue;
        else
            rise_waveform = [rise_waveform(1:idx(1)-1) ones(1,l)*amp_needed rise_waveform(idx(1):end)];
            break;
        end
    end
    
    grad_waveform = [rise_waveform rise_waveform(end:-1:1)];
else
    flat_k = (k-threshold_k);
    flat_time = flat_k/g_max/gambar;
    flat_len = floor(floor(flat_time/dt)/2)*2;
    flat_waveform = ones(1,flat_len)*g_max;
    
    flat_k = sum(flat_waveform)*dt*gambar;

    rise_time = g_max/slew_rate;
    rise_times = [0:1:ceil(rise_time/dt)]*dt;
    rise_waveform = rise_times*slew_rate;    
    rise_waveform(end)=g_max;
    rise_sum = sum(rise_waveform)*dt*gambar;
    
    k_err = (k-flat_k-rise_sum*2)/2;
    for l = 1:length(rise_waveform)
        amp_needed = k_err/(l*dt)/gambar;
        idx = find(amp_needed<rise_waveform);
        if ( length(idx) == 0 )
            continue;
        else
            rise_waveform = [rise_waveform(1:idx(1)-1) ones(1,l)*amp_needed rise_waveform(idx(1):end)];
            break;
        end
    end
    
    grad_waveform = [rise_waveform flat_waveform rise_waveform(end:-1:1)];
    
end

grad_waveform = grad_waveform(:)*sign_k;

return;