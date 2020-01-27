function [nrmse, traj_freq_kx, traj_freq_ky] = determine_pe_field_greedy3(d,sens,mask,nfreq,type,ncand, fieldmap,xfov,yfov,pulse_params,miniMax,ishard)
if ~exist('ishard')
    ishard = 0; 
end
nrmse = [];
% ncand = 5;
n_precand = ncand;

freqs.kx = zeros(nfreq,1);
freqs.ky = zeros(nfreq,1);
if ishard == 1
    durations = setup_duration_hard(freqs,xfov,yfov,pulse_params);
else
    durations = setup_duration_sinc(freqs,xfov,yfov,pulse_params);
end

r = d;
ncoils = size(sens,3);

pre_kx = [];
pre_ky = [];

freqs.kx = [];
freqs.ky = [];
pe_cands = ones(nfreq,1);

nx = size(d,2);
ny = size(d,1);
x = [ceil(-nx/2):1:floor((nx-1)/2)];
y = [ceil(-ny/2):1:floor((ny-1)/2)];
[xx,yy]=meshgrid(x,y);
d_red = d(mask);
d_norm = norm(d_red);
An = [];
x_init = [];
n_iter = ncoils*3;

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
    
    cur_kx = [];
    cur_ky = [];
    for iter = 1:ncand
        [maximum,freq_idx] = max(coherence(:));
        coherence(freq_idx)=0;
        found = 0;
        kx = floor((freq_idx-1)/ny);
        if ( kx >= nx/2 )
            kx = kx- nx;
        end
        
        ky = mod(freq_idx-1,ny);
        if ( ky >= ny/2 )
            ky = ky-ny;
        end
        
        for search = 1:length(pre_kx)
            if ( (pre_kx(search)==kx) && (pre_ky(search)==ky) )
                found = 1;
                break;
            end
        end
        
        if found == 0
            cur_kx = [cur_kx;kx];
            cur_ky = [cur_ky;ky];
        end
    end
    
    all_cands_kx = [pre_kx;cur_kx];
    all_cands_ky = [pre_ky;cur_ky];
    
    
    nrmse = zeros(size(all_cands_kx));
    x_init = [x_init;zeros(ncoils,1)];
    
    w = exp(-1i*2*pi*fieldmap*durations(sel));
    nrmse_min = 1;
    x_est_chosen = zeros(size(x_init));
    Ak_chosen = zeros(length(d_red),ncoils);
    for iter = 1:length(all_cands_kx)
        freq = exp(1i*2*pi*xx*all_cands_kx(iter)/nx+1i*2*pi*yy*all_cands_ky(iter)/ny).*w;
        Ak = [];
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
            nrmse(iter) = norm(d_red-Ank*x_est)/d_norm;
        else
            [x_est, nrmse(iter)] = minInfiniteNorm(Ank,d_red); % calculate the min Infinite Norm coefficient c.
        end
        if nrmse(iter) < nrmse_min
            x_est_chosen = x_est;
            Ak_chosen = Ak;
            if miniMax == 0
                Ank_inv_chosen = Ank_inv;
            end
            nrmse_min = nrmse(iter);
        end
    end
    
    
    [nrmse,idx] = sort(nrmse,'ascend');
    all_cands_kx = all_cands_kx(idx);
    all_cands_ky = all_cands_ky(idx);
    
    
    %{
    disp('all candidate members');
    [all_cands_kx all_cands_ky]
    keyboard;
    %%{
    disp('press enter');
    pause();
    %}
    freqs.kx = [all_cands_kx(1);freqs.kx];
    freqs.ky = [all_cands_ky(1);freqs.ky];
    pre_kx = all_cands_kx(1:n_precand);
    pre_ky = all_cands_ky(1:n_precand);
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