function [B1,gx,gy,gz,projIm] = stack_spiral_rf(d,fmap,roi,sens,xyzRange,trajPara,weightIm)
addpath ../
traj = 'spiral';
FOV = trajPara.Fov; 
FOVz = trajPara.Fovz; 
xydimRF = trajPara.ncycle; 
zcycle = trajPara.ncycleZ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a k-space trajectory and its gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load k-space trajectories
if ~isvar('k')
	%load exckspace
	gen_kspace;

	if strcmp(traj,'spiral')
	k = ksp;g = gsp;NN = NNsp;
	else
	k = kep;g = gep;NN = NNep;
	end
	
	Nt = size(k,1); % number of samples in pulses
	dt = 4e-6; % seconds, RF and gradient sampling period
end
k = flipud(k);
g = flipud(g); 


grad_info.slew_rate = 15000;  %g/cm/sec
grad_info.g_max = 4; %g/cm; 
dt = 4*10^-6; %sec
zBlip = gen_trapezoid(1/FOVz,grad_info,dt); 
zgrad_dphase = gen_trapezoid(1/FOVz*(zcycle-1)/2,grad_info,dt); 
gzero = zeros(length(zBlip),1); 
gx = g(:,1); 
gy = g(:,2); 
gz = zeros(size(g(:,1))); 
for iz = 2:zcycle
   gx = [gx; gzero; g(:,1)];
   gy = [gy; gzero; g(:,2)];
   gz = [gz; zBlip; zeros(size(g(:,1)))]; 
end
gx = [gx; zeros(size(zgrad_dphase))];
gy = [gy; zeros(size(zgrad_dphase))];
gz = [gz; -zgrad_dphase]; 


dt = 4e-3; 
gamma = 4.257;
g3d = [gx,gy,gz]; 
k3d = -flipud(cumsum(flipud(g3d))*dt*gamma); 
figure,plot3(k3d(:,1),k3d(:,2),k3d(:,3));

save k_sos g3d k3d
% [B1,gx,gy,gz,projIm] = cont_gen_rf(d,fmap,roi,sens,xyzRange,k3d, g3d); 
[B1,gx,gy,gz,projIm] = cont_gen_rf_exactA(d,fmap,roi,sens,xyzRange, weightIm, k3d); 

