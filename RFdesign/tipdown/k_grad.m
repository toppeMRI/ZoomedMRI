% calculate gradient over kx, ky, kz
function [dkx, dky, dkz, r] = k_grad(k, d, sens, mask,weightIm)
nx = size(d,1);
ny = size(d,2);
nz = size(d,3);
x = [ceil(-nx/2):1:floor((nx-1)/2)];
y = [ceil(-ny/2):1:floor((ny-1)/2)];
z = [ceil(-nz/2):1:floor((nz-1)/2)];
[xx,yy,zz]=ndgrid(x,y,z);

kx = k(:,1); 
ky = k(:,2); 
kz = k(:,3); 

ncoils = size(sens,4); 
%% construct the matrix A
% w = exp(-1i*2*pi*fieldmap*durations(sel));

A = []; 
for ik = 1:length(kx)
   freq = exp(1i*2*pi*xx*kx(ik)/nx+1i*2*pi*yy*ky(ik)/ny...
      +1i*2*pi*zz*kz(ik)/nz); % .*w;
   Ak = [];
   for ic = 1:ncoils
      fic = freq.*squeeze(sens(:,:,:,ic)).*sqrt(weightIm);
      Ak = [Ak fic(mask)];
   end
   A = [A Ak];
end

b = pinv(A)*d(mask); 
% calculate error
% a = embed(A*b, mask); 
r = d - embed(A*b, mask); 

%% Compute the gradient
Db = diag_sp(b');
X = diag_sp(xx(mask));
Y = diag_sp(yy(mask)); 
Z = diag_sp(zz(mask)); 
Sr_tmp = sens(:,:,:,1); 
Sr = diag_sp(Sr_tmp(mask)');
W = diag_sp(weightIm(mask));
e = r(mask); 

dkx = 2*real(1i*Db*A'*X*Sr*W*e); 
dky = 2*real(1i*Db*A'*Y*Sr*W*e); 
dkz = 2*real(1i*Db*A'*Z*Sr*W*e); 

