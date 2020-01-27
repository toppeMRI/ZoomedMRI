function [nrmse, traj_kx, traj_ky, traj_kz] = determine_pe_field_inverseA_3d(d,mask,nfreq)

% x = [ceil(-nx/2):1:floor((nx-1)/2)];
% y = [ceil(-ny/2):1:floor((ny-1)/2)];
% z = [ceil(-nz/2):1:floor((nz-1)/2)];
% [xx,yy,zz]=ndgrid(x,y,z);
% kxCand = [ceil(-nx/2):1:floor((nx-1)/2)];
% kyCand = [ceil(-ny/2):1:floor((ny-1)/2)];
% kzCand = [ceil(-nz/2):1:floor((nz-1)/2)];
% [xxk, kxx] = ndgrid(xx(:),kxCand(:));
% [yyk, kyy] = ndgrid(yy(:),kyCand(:));
% [zzk, kzz] = ndgrid(zz(:),kzCand(:));
% 
[nx,ny,nz] = size(mask); 
A = Gdft('mask',mask,'fftshift',0); 
wts_dft = A*d; 
% wts_embed = embed(wts_dft,true(size(mask))); 
% figure,im(wts_embed,'.');  
[wts_sort, sort_ind] = sort(wts_dft(:),'descend');
sort_first = sort_ind(1:300);
wts_thre = zeros(size(wts_dft)); 
wts_thre(sort_first) = 1; 
figure,im(wts_thre); 

for ifreq = 1:nfreq
    freq_idx = sort_ind(ifreq);
    
    kz = floor((freq_idx-1)/nx/ny);
    kz_tmp = kz;
    if ( kz >= nz/2 )
        kz = kz- nz;
    end
    traj_kz(ifreq) = kz; 
    
    ky = floor((freq_idx-1-kz_tmp*nx*ny)/ny);
    if ( ky >= ny/2 )
        ky = ky-ny;
    end
    traj_ky(ifreq) = ky; 
    
    kx = mod((freq_idx-1-kz_tmp*nx*ny),nx);
    if ( kx >= nx/2 )
        kx = kx-nx;
    end
    traj_kx(ifreq) = kx; 
end
nrmse = 1; 
traj_kx = flipud(traj_kx'); 
traj_ky = flipud(traj_ky'); 
traj_kz = flipud(traj_kz'); 
end