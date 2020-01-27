function [nrmse, traj_freq_kx, traj_freq_ky] = determine_pe_field_greedy3_hao(d,sens,mask,nfreq,type,ncand, fieldmap,xfov,yfov,pulse_params,miniMax)

nrmse = [];
% ncand = 5;

freqs.kx = zeros(nfreq,1);
freqs.ky = zeros(nfreq,1);
durations = setup_duration(freqs,xfov,yfov,pulse_params);


r = d;
ncoils = size(sens,3);


freqs.kx = [];
freqs.ky = [];

nx = size(d,2);
ny = size(d,1);
x = [ceil(-nx/2):1:floor((nx-1)/2)];
y = [ceil(-ny/2):1:floor((ny-1)/2)];
[xx,yy]=meshgrid(x,y);
d_red = d(mask);
d_norm = norm(d_red);
An = [];
x_init = [];

for sel = nfreq:-1:1
    r_n = r.*exp(1i*2*pi*fieldmap*durations(sel));
    coherence = compute_coherence(r_n, sens, size(sens,3), type);
    %%{
    [sorted, sorted_idx] = sort(coherence(:),'descend');
    sorted_ipt = coherence(sorted_idx);
    kx_sorted = zeros(length(sorted_idx),1);
    ky_sorted = zeros(length(sorted_idx),1);
    for ll = 1:length(sorted_idx)
        kx_sorted(ll) = floor((sorted_idx(ll)-1)/ny);
        if ( kx_sorted(ll) >= nx/2 )
            kx_sorted(ll) = kx_sorted(ll)-nx;
        end
        
        ky_sorted(ll) = mod(sorted_idx(ll)-1,ny);
        if (ky_sorted(ll) >= ny/2 )
            ky_sorted(ll) = ky_sorted(ll)-ny;
        end
    end
    %disp('new candidate members');
    %[kx_sorted(1:10) ky_sorted(1:10)]
    %%}
    
    
    
    
    nrmse = zeros(ncand);
    x_init = [x_init;zeros(ncoils,1)];
    
    w = exp(-1i*2*pi*fieldmap*durations(sel));
    nrmse_min = 1;
    x_est_chosen = zeros(size(x_init));
    Ak_chosen = zeros(length(d_red),ncoils);
    
    for icand = 1:ncand
        Ak = []; 
        freq = exp(1i*2*pi*xx*kx_sorted(icand)/nx+1i*2*pi*yy*ky_sorted(icand)/ny).*w;
       
        for ic = 1:ncoils
            fic = freq.*squeeze(sens(:,:,ic));
            Ak = [Ak fic(mask)];
        end
        Ank = [An Ak];
        
        if miniMax == 0
            if sel == nfreq
                Ank_inv = inv(Ak'*Ak);
            else
                uk = An'*Ak;
                Qk = inv(Ak'*Ak-uk'*An_inv*uk);
                Ank_inv = [An_inv+An_inv*uk*Qk*uk'*An_inv -An_inv*uk*Qk;-Qk*uk'*An_inv Qk];
            end
            
            x_est = (Ank_inv*(Ank'*d_red));
            nrmse(icand) = norm(d_red-Ank*x_est)/d_norm;
        else
            [x_est, nrmse(icand)] = minInfiniteNorm(Ank,d_red); % calculate the min Infinite Norm coefficient c.
        end
        if nrmse(icand) < nrmse_min
            x_est_chosen = x_est;
            Ak_chosen = Ak;
            if miniMax == 0
                Ank_inv_chosen = Ank_inv;
            end
            nrmse_min = nrmse(icand);
        end
    end
    
    
    [nrmse,idx] = sort(nrmse,'ascend');
    
    
    %{
    disp('all candidate members');
    [all_cands_kx all_cands_ky]
    keyboard;
    %%{
    disp('press enter');
    pause();
    %}
    freqs.kx = [kx_sorted(idx(1));freqs.kx];
    freqs.ky = [ky_sorted(idx(1));freqs.ky];

    %keyboard;
    x_init = x_est_chosen;
    An = [An Ak_chosen];
    if miniMax == 0
        An_inv = Ank_inv_chosen;
    end
    r = d - embed(An*x_init,mask);
end

traj_freq_kx = freqs.kx;
traj_freq_ky = freqs.ky;
nrmse = nrmse(1);

end