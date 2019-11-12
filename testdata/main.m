
% from Ray
load P51200_res.mat

% zero-pad to isotropic matrix
P = sqrt(sum(abs(fPS(:,:,:,:,1)).^2,4));
kP = fftshift(fftn(fftshift(P)));
kP2 = zeros(80,80,80);
kP2(:,:,(end/2-20):(end/2+19)) = kP;
P2 = fftshift(ifftn(fftshift(kP2)));
P2 = abs(P2);

% create Localizer.h5 that can be loaded into ~/github/toppeMRI/SlicePlanner/GUI/
addpath ~/github/toppeMRI/SlicePlanner/LocalizerScan/   % im2hdf
%im2hdr(P2); 

% object support
roi = getROI('ROI.h5',1);
mask = roi2mask(roi,80,4);
objSupportMask = mask & (P2>0.15*max(P2(:)));
save objSupportMask objSupportMask

% IV mask
roi = getROI('ROI.h5',2);
iv = roi2mask(roi,80,4);
save iv iv
