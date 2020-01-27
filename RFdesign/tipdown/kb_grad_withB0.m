function [dkx, dky, dkz, db, e, As, J] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt, xyzRange)
% calculate gradient over kx, ky, kz, and b
% function [dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt, xyzRange)

nx = size(d,1);
ny = size(d,2);
nz = size(d,3);
if ~exist('xyzRange')
   x = [ceil(-nx/2):1:floor((nx-1)/2)];
   y = [ceil(-ny/2):1:floor((ny-1)/2)];
   z = [ceil(-nz/2):1:floor((nz-1)/2)];
   defaultXYZ = 1; 
else
   x = xyzRange.x/(xyzRange.x(2)-xyzRange.x(1)); 
   y = xyzRange.y/(xyzRange.y(2)-xyzRange.y(1)); 
   z = xyzRange.z/(xyzRange.z(2)-xyzRange.z(1)); 
   defaultXYZ = 0; 
end
[xx,yy,zz]=ndgrid(x,y,z);

kx = k(:,1); 
ky = k(:,2); 
kz = k(:,3); 

ncoils = size(sens,4); 
%% construct the matrix A and evaluate residule
if defaultXYZ == 1
   A = formA(k, sens, mask, fieldmap, tt); 
else
   A = formA(k, sens, mask, fieldmap, tt, xyzRange); 
end

% combine sens and A into As
if ncoils ==1
   As = A;
else
   Atmp = [];
   for ic = 1:ncoils
      sensi = sens(:,:,:,ic);
      Atmp = [Atmp; spdiag(sensi(mask))*A];
   end
   As = reshape(Atmp, [sum(mask(:)),ncoils*length(kx)]);
end
    
if b == 0  % calculate an initialization for b if set b = 0 as the input . 
   b = dpls(As,d,beta,mask);
end

% calculate error
% a = embed(A*b, mask); 
% r = d - embed(A*b, mask); 
% e = r(mask); 
e = d(mask) - As*b;

%% Compute the gradients and Jacobian
Db = diag_sp(b); % has transpose here
X = diag_sp(xx(mask));
Y = diag_sp(yy(mask)); 
Z = diag_sp(zz(mask)); 
Sr_tmp = sens(:,:,:,1); 
Sr = diag_sp(Sr_tmp(mask)); 
W = diag_sp(weightIm(mask));

Jx = Sr*X*A*Db; % remove the negative sign here
Jy = Sr*Y*A*Db; 
Jz = Sr*Z*A*Db; 
J = [Jx, Jy, Jz]; 

dkx = 2*real(1i*Jx'*W*e)/nx; 
dky = 2*real(1i*Jy'*W*e)/ny; 
dkz = 2*real(1i*Jz'*W*e)/nz; 

db = A'*(A*b-d(mask)); 

