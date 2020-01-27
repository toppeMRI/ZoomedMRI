% Simulate reduced FOV acquisition
% To generate the spiral k-space traj:
% jfnielse@brooks:~/projects/psd/umpsd/matlab/spiral/vds.m

%% load input image and k-space
close all; 
clear; 
% function imRFOV = simuRFOV(inputIm)
% inputIm = d;
% inputIm = M_bloch/sin(10/180*pi); 

%%%%%%%%%% for rFOV paper
load rFOV_data_oct23_scanJon
inputIm = imbssfp; 
inputImtmp = []; 
for iz = 1:size(inputIm,3)
inputImtmp(:,:,iz) = imresize(inputIm(:,:,iz),[240,240]); 
end
inputIm = inputImtmp; 


%%%%%%%%%% for images with cranial nerve 
% load rfov_cranialNerves_fromJon
% inputIm = im_iac_iv; 
% centerR = 86; 
% centerC = 133;
% inputIm = circshift(inputIm, [120-centerR, 120-centerC]); 

[nx, ny, nz] = size(inputIm);
fov1 = 24; % fov of object 
fov2 = 8; % fov of recon
% load k_downsample_spiral_jon; % 10 cm fov, 33 nleafs, k in cycle/fov
% kx = kx/fov2; 
% ky = ky/fov2; 
load rFOV_ktraj_24-12-8cm.mat; % 24/12/8 cm fov, 36/18/12 nleafs, k in cycle/cm
k = genvarname(['k' num2str(fov2) 'all']); 
% k = ['k' num2str(fov2) 'all']; 
eval(['kx = real(' k ');']); 
eval(['ky = imag(' k ');']); 
kxAll = kx(:); 
kyAll = ky(:); 
k = [kxAll, kyAll]; 

nt = size(k,1); 
dt = 4e-6; 
gambar = 4257; %Hz/G
gam = gambar*2*pi;

%% from the magnetization to k-space
% A_read = formA_readout(k,nx,ny,fov1);  % nufft in kx-ky
A_read = gam*dt*Gmri(k, true(nx,ny), 'fov', [fov1; fov1]); 
% fov1 = fov1*nz/(nz-1); 
% nufft_args = {[nx ny],[6 6],[nx ny]*2,[nx-1 ny-1]/2,'minmax:kb'};
% A_read = gam*dt*Gmri(k, true(nx,ny), 'fov', [fov1; fov1],'nufft',nufft_args); 

for iz = 1:nz
    tmp = inputIm(:,:,iz); 
    inputIm_nufft(:,iz) = A_read*tmp(:);  
end
kdata = fft(inputIm_nufft,[],2);  % fft in z

%% undersample in z (set every other kz slice to zero)
downFactor = 6; 
kzslices = 1:downFactor:nz; 
kdataUnder = zeros(size(kdata));
kdataUnder(:,kzslices) = kdata(:,kzslices); 

%% recon
reso = fov1/nx; 
nx2 = round(fov2/reso); 
ny2 = round(fov2/reso);  
fov2 = nx2*reso; % resolution have to be matched 
nz2 = round(nz/downFactor); 
% A_read2 = formA_readout(k,nx2,ny2,fov2);
A_read2 = gam*dt*Gmri(k, true(nx2,ny2), 'fov', [fov2; fov2]); 
W = diag_sp(ones(length(k),1)); 
kdata_inufft = zeros(nx2,ny2,nz); 
for iz = 1:length(kzslices)
    tmp = qpwls_pcg1_hao(zeros(nx2*ny2,1), A_read2, W, kdataUnder(:,kzslices(iz)), 0, 'niter', 20); 
    kdata_inufft(:,:,kzslices(iz)) = reshape(tmp, [nx2, ny2]); 
end
outputAll = ifft(kdata_inufft,[],3)*downFactor; 
zslices = floor(nz/2)-floor(nz2/2)+1:floor(nz/2)+ceil(nz2/2); 
output = outputAll(:,:,zslices); 

%% plot
nz2 = length(zslices); 
zr = 1:nz2; 
zr = floor(length(zslices)/2)-1:floor(length(zslices)/2)+4; % for 6 slices 
figure,im(output(:,:,zr)),%colormap default; 
title(sprintf('rFOV acquisition; down sample fold: %.1f (x-y), %d (z)',fov1/fov2, downFactor)); 
colorbar; 
% print('-depsc', sprintf('rFOV_%.1f_x-y_%d_z.eps',fov1/fov2,downFactor)); 

cropIm = inputIm(nx/2-floor(nx2/2)+1:nx/2+ceil(nx2/2), ny/2-floor(ny2/2)+1:ny/2+ceil(ny2/2),zslices);
figure,im(cropIm(:,:,zr)),%colormap default;
title('cropped image'); 
colorbar; 
% print('-depsc',  sprintf('cropIm_%d_xyFOV_%d_z.eps',fov2,downFactor)); 
% 
diffIm = output - cropIm; 
figure, im(diffIm(:,:,zr)),% colormap default; 
colorbar; 
% print('-depsc',  sprintf('diffIm_%d_xyFOV_%d_z.eps',fov2,downFactor)); 

figure, im('row',3,cat(3,cropIm(:,:,zr),output(:,:,zr),5*diffIm(:,:,zr)),[0 max(abs(cropIm(:)))]); 
colorbar; 
title(sprintf('FOVxy: %d cm; FOVz: %.1f cm',fov2, 21/downFactor));
% print('-depsc', sprintf('allIm_%d_xyFOV_%d_z.eps',fov2,downFactor)); 

%% calculate error vs downsampling
xr = nx2/2-6/fov2*nx2/2 : nx2/2+6/fov2*nx2/2-1;
yr = xr; 
roi = false(nx2,ny2,nz2); 
roi(xr,yr,zr) = true; 
figure,im(roi); 
diffIm_norm = norm(diffIm(roi)); 
cropIm_norm = norm(cropIm(roi)); 
output_norm = norm(output(roi)); 

rmse = diffIm_norm./cropIm_norm 
save(sprintf('rFOV_bssfp_%d_x-y_%d_z',fov2,downFactor), 'diffIm_norm',...
   'cropIm_norm', 'output_norm','rmse'); 
% end
