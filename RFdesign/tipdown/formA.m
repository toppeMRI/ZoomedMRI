function [A,mask,xx,yy,zz,nx,ny,nz] = formA(k, sens, mask,fieldmap,tt, xyzRange)
% construct the system matrix
% [A] = formA(k, sens, mask,fieldmap,tt, xyzRange)
% form x, y, z if not as input
% k is in the unit of cycle/fov
% x,y,z are in grid, not cm
% tt in sec

nx = size(fieldmap,1);
ny = size(fieldmap,2);
nz = size(fieldmap,3);

if ~exist('xyzRange')
   x = [ceil(-nx/2):1:floor((nx-1)/2)];
   y = [ceil(-ny/2):1:floor((ny-1)/2)];
   z = [ceil(-nz/2):1:floor((nz-1)/2)];
else
   x = xyzRange.x/(xyzRange.x(2)-xyzRange.x(1));
   y = xyzRange.y/(xyzRange.y(2)-xyzRange.y(1));
   z = xyzRange.z/(xyzRange.z(2)-xyzRange.z(1));
end
[xx,yy,zz]=ndgrid(x,y,z);

kx = k(:,1);
ky = k(:,2);
kz = k(:,3);

ncoils = size(sens,4);

%% construct the matrix A and evaluate residule

% w = exp(-1i*2*pi*fieldmap*durations(sel));
% A = zeros(sum(mask(:)),ncoils*length(kx));
% for ik = 1:length(kx)
%    freq = exp(1i*2*pi*xx*kx(ik)/nx+1i*2*pi*yy*ky(ik)/ny...
%       +1i*2*pi*zz*kz(ik)/nz); % .*w;
%    % Ak = zeros(sum(mask(:)),ncoils);
%    for ic = 1:ncoils
%       fic = freq.*squeeze(sens(:,:,:,ic));
%       % Ak(:,ic) = fic(mask);
%       A(:,(ik-1)*ncoils+ic) = fic(mask);
%    end
%    % A(:,ncoils*(ik-1)+1:ncoils*ik) = Ak;
% end

gambar = 4257; %Hz/G
gam = gambar*2*pi;
dt = 4e-6; % second

Ab0 = 1i*2*pi*fieldmap(mask(:))*tt;
%freq = 1i*gam*dt*exp(Ab0 + 1i*2*pi*(xx(mask(:))*(kx'/nx) + yy(mask(:))*(ky'/ny) + zz(mask(:))*(kz'/nz)));
A = 1i*gam*dt*exp(Ab0 + 1i*2*pi*(xx(mask(:))*(kx'/nx) + yy(mask(:))*(ky'/ny) + zz(mask(:))*(kz'/nz)));

% don't put the sensitivity matrix in to system matrix A
%A = reshape(freq, [sum(mask(:)),length(kx)]);

% if ncoils ==1
%    As = A;
% else
%    Atmp = [];
%    for ic = 1:ncoils
%       sensi = sens(:,:,:,ic);
%       Atmp = [Atmp; spdiag(sensi(mask))*freq];
%    end
%    As = reshape(Atmp, [sum(mask(:)),ncoils*length(kx)]);
% end
