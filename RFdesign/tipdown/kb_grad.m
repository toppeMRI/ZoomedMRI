% calculate gradient over kx, ky, kz, and b
function [dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm,beta)
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

freq = exp(1i*2*pi*xx(mask)*kx'/nx + 1i*2*pi*yy(mask)*ky'/ny + 1i*2*pi*zz(mask)*kz'/nz);
Atmp = [];
for ic = 1:ncoils
   sensi = sens(:,:,:,ic);
   Atmp = [Atmp; spdiag(sensi(mask))*freq];
end
A = reshape(Atmp, [sum(mask(:)),ncoils*length(kx)]); 


if b == 0  % calculate an initialization for b if set b = 0 as the input . 
   b = dpls(A,d,beta,mask);
end
% calculate error
% a = embed(A*b, mask); 
% r = d - embed(A*b, mask); 
% e = r(mask); 
e = d(mask) - A*b;
%% Compute the gradients
Db = diag_sp(b'); % has transpose here
X = diag_sp(xx(mask));
Y = diag_sp(yy(mask)); 
Z = diag_sp(zz(mask)); 
Sr_tmp = sens(:,:,:,1); 
Sr = diag_sp(Sr_tmp(mask)'); % has transpose here
W = diag_sp(weightIm(mask));


dkx = 2*real(1i*Db*A'*X*Sr*W*e); 
dky = 2*real(1i*Db*A'*Y*Sr*W*e); 
dkz = 2*real(1i*Db*A'*Z*Sr*W*e); 

db = A'*(A*b-d(mask)); 

