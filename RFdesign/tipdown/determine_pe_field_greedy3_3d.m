function [nrmse, traj_freq_kx, traj_freq_ky, traj_freq_kz] = determine_pe_field_greedy3_3d(d,sens,mask,nfreq,type,ncand,fieldmap,xfov,yfov,zfov,pulse_params,miniMax,ishard,isCont)
% modified from daehyun's code : determine_pe_field_greedy
d(~mask) = 0; % remove the effect of outside ROI in calculation correlation between residual and system matrix
if ~exist('ishard','var')
   ishard = 0;
end

if ~exist('isCont','var')
   isCont = 1; % for continue traj design. Ignore the maximum phase jump constrain.
   % doesn't allow the same phase encoding locations be selected again
end

load weighting; 
nrmse = [];
% ncand = 5;
n_precand = ncand;

freqs.kx = zeros(nfreq,1);
freqs.ky = zeros(nfreq,1);
freqs.kz = zeros(nfreq,1);
if ishard == 1
   durations = setup_duration_hard(freqs,xfov,yfov,zfov,pulse_params);
else
   durations = setup_duration_sinc(freqs,xfov,yfov,pulse_params);
end


r = d;
ncoils = size(sens,4);

pre_kx = [];
pre_ky = [];
pre_kz = [];

freqs.kx = [];
freqs.ky = [];
freqs.kz = [];
pe_cands = ones(nfreq,1);

nx = size(d,1);
ny = size(d,2);
nz = size(d,3);

x = [ceil(-nx/2):1:floor((nx-1)/2)];
y = [ceil(-ny/2):1:floor((ny-1)/2)];
z = [ceil(-nz/2):1:floor((nz-1)/2)];
[xx,yy,zz]=ndgrid(x,y,z);
d = d.*sqrt(weightIm); % weighted least square; 
d_red = d(mask);
d_norm = norm(d_red);
An = [];
x_init = [];
n_iter = ncoils*3;

for sel = nfreq:-1:1
   r_n = r.*exp(1i*2*pi*fieldmap*durations(sel));
   coherence = compute_coherence_3d(r_n, sens, size(sens,4), type);
   if sel == nfreq - 12
      sel;
   end 
   %%{
   [sorted, sorted_idx] = sort(coherence(:),'descend');
   sorted_ipt = coherence(sorted_idx);
   %     kx_sorted = zeros(length(sorted_idx),1);
   %     ky_sorted = zeros(length(sorted_idx),1);
   %     kz_sorted = zeros(length(sorted_idx),1);
   %     for ll = 1:length(sorted_idx)
   %         kz_sorted(ll) = floor((sorted_idx(ll)-1)/nx/ny);
   %         kz_sorted_tmp = kz_sorted(ll);
   %         if ( kz_sorted(ll) >= nz/2 )
   %             kz_sorted(ll) = kz_sorted(ll)-nz;
   %         end
   %
   %         ky_sorted(ll) = floor((sorted_idx(ll)-1-nx*ny*kz_sorted_tmp)/nx);
   %         if (ky_sorted(ll) >= ny/2 )
   %             ky_sorted(ll) = ky_sorted(ll)-ny;
   %         end
   %
   %         kx_sorted(ll) = mod((sorted_idx(ll)-1-nx*ny*kz_sorted_tmp),nx);
   %         if (kx_sorted(ll) >= nx/2 )
   %             kx_sorted(ll) = kx_sorted(ll)-nx;
   %         end
   %     end
   %disp('new candidate members');
   %[kx_sorted(1:10) ky_sorted(1:10)]
   %%}
   
   cur_kx = [];
   cur_ky = [];
   cur_kz = [];
   for iter = 1:ncand
      [maximum,freq_idx] = max(coherence(:));
      coherence(freq_idx)=0; % exclude the selected PE locations
      found = 0;
      kz = floor((freq_idx-1)/nx/ny);
      kz_tmp = kz;
      if ( kz >= nz/2 )
         kz = kz- nz;
      end
      
      ky = floor((freq_idx-1-kz_tmp*nx*ny)/ny);
      if (ky> 100)
         ky
         freq_idx
      end
      if ( ky >= ny/2 )
         ky = ky-ny;
      end
      
      kx = mod((freq_idx-1-kz_tmp*nx*ny),nx);
      if ( kx >= nx/2 )
         kx = kx-nx;
      end
      
      for search = 1:length(pre_kx)
         if ( (pre_kx(search)==kx) && (pre_ky(search)==ky) && (pre_kz(search)==kz))
            found = 1;
            break;
         end
      end
      
      if found == 0
         cur_kx = [cur_kx;kx]; % new kx selected from inner product is not appeared in previous top cands
         cur_ky = [cur_ky;ky];
         cur_kz = [cur_kz;kz];
      end
   end % end of selecting ncand cands, those cands are combined with previous cands to form a new cands set
   
   all_cands_kx = [pre_kx;cur_kx];
   all_cands_ky = [pre_ky;cur_ky];
   all_cands_kz = [pre_kz;cur_kz];
   
   if isCont == 1
      if sel == nfreq - 1 % In the second run of sel loop, the k location corresponding to An (0 0 0) are in pre_kx, and thus all_cands,
         % which leads to singular matrix in the following loop (e.g. Ank = [A(0) A(0)]).
         % This only appear once and should not cause problem even
         % if I didn't exclude such a k location from cands.
         all_cands_kx(1) = [];
         all_cands_ky(1) = [];
         all_cands_kz(1) = [];
      end
   end
   
   nrmse = zeros(size(all_cands_kx));
   x_init = [x_init;zeros(ncoils,1)];
   
   w = exp(-1i*2*pi*fieldmap*durations(sel));
   nrmse_min = 2;
   x_est_chosen = zeros(size(x_init));
   Ak_chosen = zeros(length(d_red),ncoils);
   for iter = 1:length(all_cands_kx)
      freq = exp(1i*2*pi*xx*all_cands_kx(iter)/nx+1i*2*pi*yy*all_cands_ky(iter)/ny...
         +1i*2*pi*zz*all_cands_kz(iter)/nz).*w;
      Ak = [];
      for ic = 1:ncoils
         fic = freq.*squeeze(sens(:,:,:,ic)).*sqrt(weightIm);
         Ak = [Ak fic(mask)];
      end
      Ank = [An Ak];
      
      if miniMax == 0
         if sel == nfreq || sel == nfreq - 1 % nfreq-1 is added by Hao
            %               Ank_inv = inv(Ak'*Ak);
            Ank_inv = inv(Ank'*Ank);
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
      
      
      if length(freqs.kx) >= 1   % Hao: only select the candidate with PE location jump smaller than threshold.
         kx_jump = (all_cands_kx(iter) - freqs.kx(1))/xfov;
         ky_jump = (all_cands_ky(iter) - freqs.ky(1))/yfov;
         kz_jump = (all_cands_kz(iter) - freqs.kz(1))/zfov;
      else
         kx_jump = 0;
         ky_jump = 0;
         kz_jump = 0;
      end
      kxyz_jump_max = pulse_params.max_pe_jump;
      if isCont == 1
         kxyz_jump_max = kxyz_jump_max*100;
      end
      jumpSmall = 0; 
      if abs(kx_jump) < kxyz_jump_max...
            && abs(ky_jump) < kxyz_jump_max...
            && abs(kz_jump) < kxyz_jump_max
         jumpSmall = 1; 
      end
      
      if nrmse(iter) < nrmse_min && jumpSmall == 1 % Hao: the selected Ank
         x_est_chosen = x_est;
         Ak_chosen = Ak;
         if miniMax == 0
            Ank_inv_chosen = Ank_inv;
         end
         nrmse_min = nrmse(iter);
         k_ind_chosen = iter; 
      end
   end % end of evaluating all cands to minimize ||[An, cands(:,iter)]x - r||
   
   freqs.kx = [all_cands_kx(k_ind_chosen);freqs.kx]; % update selected k location 
   freqs.ky = [all_cands_ky(k_ind_chosen);freqs.ky];
   freqs.kz = [all_cands_kz(k_ind_chosen);freqs.kz];
   
   x_init = x_est_chosen; % update A, x, and residual r
   An = [An Ak_chosen];
   if miniMax == 0
      An_inv = Ank_inv_chosen;
   end
   r = d - embed(An*x_init,mask);
   
   [nrmse,idx] = sort(nrmse,'ascend'); % keep some cands for next iteration
   all_cands_kx = all_cands_kx(idx);
   all_cands_ky = all_cands_ky(idx);
   all_cands_kz = all_cands_kz(idx); 
   pre_kx = all_cands_kx(1:n_precand);
   pre_ky = all_cands_ky(1:n_precand);
   pre_kz = all_cands_kz(1:n_precand);
end


traj_freq_kx = freqs.kx;
traj_freq_ky = freqs.ky;
traj_freq_kz = freqs.kz;
nrmse = nrmse(1);

%% gradient descent
% k = [freqs.kx,freqs.ky,freqs.kz]; 
% [k, cost] = GD_ktraj(k, d, sens, mask, weightIm); 

end