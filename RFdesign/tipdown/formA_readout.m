function [A] = formA_readout(k,nx,ny,fov)
% construct the system matrix: from magnetization to k-space
% [A] = formA_readout(k, nx, ny)
% k: cycle/cm

x = [ceil(-nx/2):1:floor((nx-1)/2)]/nx*fov;
y = [ceil(-ny/2):1:floor((ny-1)/2)]/ny*fov;
[xx,yy] = ndgrid(x,y);

kx = k(:,1);
ky = k(:,2);

%% construct the matrix A and evaluate residule



gambar = 4257; %Hz/G
gam = gambar*2*pi;
dt = 4e-6; % second

% Ab0 = -1i*2*pi*fieldmap(mask)*tt;
freq = 1i*gam*dt*exp(-1i*2*pi*((kx)*xx(:)' + (ky)*yy(:)'));

% don't put the sensitivity matrix in to system matrix A
A = reshape(freq, [length(kx),nx*ny]);
