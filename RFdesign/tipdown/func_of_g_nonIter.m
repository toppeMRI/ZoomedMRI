function [fval] = func_of_g_nonIter(g_coeff, d,fmap,roi,sens,xyzRange,Trf)
% [B1,gx,gy,gz,projIm] = func_of_g(g_coeff, d,fmap,roi,sens,xyzRange,Trf)
% Based on: (1) Spiral RF pulse design, small-tip (Grissom MRM2005)
%           (2) Generate RF waveforms given k-traj
%
% Units: 
% k-traj: cycle/cm
% gradient: Gauss/cm
% B1: Gauss
%
% Hao: Jun.25. 2012 Part of the code is borrow from spiralrf.m (Grissom)
%close all; 

beta = 8; % regularizar on b

Nt = round(Trf/4e-6); 
tfun = 1:Nt; tfun = tfun';
ncycle = length(g_coeff)/3; % degree of freedom
    
if 1 % fourier coeff.
    icycle = 0:ncycle-1;
    contGx = cos(2*pi*(tfun*icycle)/Nt)*g_coeff(1:ncycle);
    contGy = cos(2*pi*(tfun*icycle)/Nt)*g_coeff(ncycle+1:2*ncycle);
    contGz = cos(2*pi*(tfun*icycle)/Nt)*g_coeff(2*ncycle+1:3*ncycle);
else
    contGx = polyval(g_coeff(1:ncycle),tfun);
    contGy = polyval(g_coeff(ncycle+1:2*ncycle),tfun);
    contGz = polyval(g_coeff(2*ncycle+1:3*ncycle),tfun);
end
contG = [contGx, contGy, contGz];

dt = 4e-3; 
gamma = 4.257; 
k = -flipud(cumsum(flipud(contG),1)*dt*gamma); 
figure(11),plot(contG); 
%tipangle = 16; % pulse flip angle
traj = 'spins';
if ~isvar('reverseK')
    reverseK = 0; 
end
Nz = size(fmap,3); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('sens')
	%load fdtdsens
	if 0
		load sens;  % from ~/projects/parallelrf/data/Apr2011/9/8coil_sens.mat
		%mask = sum(abs(sens),3) > 0.1*max(col(sum(abs(sens),3))) & roi;
	else
		sens = ones(1,64,64,Nz);
    end
end

% Hao: reshape the coil dimension as the 4th dimension
sens_tmp = [];
if size(sens,1) < 9
   for ic = 1:size(sens,1)
      sens_tmp(:,:,:,ic) = sens(ic,:,:,:);
   end
elseif length(size(sens)) < 4
   sens_tmp = ones(64,64,Nz,1); 
end
sens = sens_tmp; 

%% load mask; 
%mask = logical(ones(size(roi)));
%sens(:,:,2:2:end) = [];
%save mask mask
load mask; 
load roi; 
dim = length(xyzRange.x); % dimension of square x-y grid
dimz = length(xyzRange.z); 
Nc = size(sens,4); % number of Tx coils
FOV = -xyzRange.x(1)*2; 
FOVz = -xyzRange.z(1)*2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the desired flip angle pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xx, yy, zz] = ndgrid(xyzRange.x,xyzRange.y,xyzRange.z); 
sens = sens/max(abs(sens(:)));


k_int(:,1) = k(:,1)*FOV; 
k_int(:,2) = k(:,2)*FOV; 
k_int(:,3) = k(:,3)*FOVz; 
k = k_int; 

[dkx, dky, dkz, db, e, A] = kb_grad(k, 0, d, sens, mask,ones(size(mask)),beta);
b = dpls(A,d,beta,mask); 
projIm = embed(A*b,mask); 
M_diff = projIm - d;
fval = norm(M_diff(mask)); 

return;

