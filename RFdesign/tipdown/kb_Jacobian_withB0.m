function [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt, xyzRange)
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

Nt = size(k,1);
ncoils = size(sens,4); 
Np = length(xx(mask));
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

% Ab0 = 1i*2*pi*fieldmap(mask)*tt;
% freq = exp(Ab0 + 1i*2*pi*xx(mask)*kx'/nx + 1i*2*pi*yy(mask)*ky'/ny + 1i*2*pi*zz(mask)*kz'/nz);
% Atmp = [];
% for ic = 1:ncoils
%    sensi = sens(:,:,:,ic);
%    Atmp = [Atmp; spdiag(sensi(mask))*freq];
% end
% A = reshape(Atmp, [sum(mask(:)),ncoils*length(kx)]);
if defaultXYZ == 1
   A = formA(k, sens, mask, fieldmap, tt);
else
   A = formA(k, sens, mask, fieldmap, tt, xyzRange);
end

% combine sens and A into As
Atmp = [];
for ic = 1:ncoils
   sensi = sens(:,:,:,ic);   
   sensiv = sensi(mask); 
   sensMat = repmat(sensiv,[1,Nt]);
   Atmp = [Atmp; spdiag(sensi(mask))*A];
%    Atmp = [Atmp; sensMat.*A];
end
As = reshape(Atmp, [sum(mask(:)),ncoils*length(kx)]);


if b == 0  % calculate an initialization for b if set b = 0 as the input .
   b = dpls(As,d,beta,mask);
end

% calculate error
% a = embed(A*b, mask);
% r = d - embed(A*b, mask);
% e = r(mask);
e = d(mask) - As*b;

%% Compute the gradients and Jacobian


xv = xx(mask);
yv = yy(mask);
zv = zz(mask);

X = spdiag(xx(mask),0,Np,Np);
Y = spdiag(yy(mask),0,Np,Np);
Z = spdiag(zz(mask),0,Np,Np);
% W = spdiag(weightIm(mask),0,Np,Np);

Jx = zeros(Np, Nt);
Jy = zeros(Np, Nt);
Jz = zeros(Np, Nt);

for ic = 1:ncoils
   Db = spdiag(b((ncoils-1)*Nt+1:(ncoils)*Nt),0,Nt,Nt);
   Sr_tmp = sens(:,:,:,ic);
   
   % % sparse matrix only works with double
   %    Sr = spdiag(Sr_tmp(mask),0,Np,Np);
   %    Jx = Jx + Sr*X*double(A)*Db; % remove the negative sign here
   %    Jy = Jy + Sr*Y*double(A)*Db;
   %    Jz = Jz + Sr*Z*double(A)*Db;
   
   % % avoid sparse matrix because it is not compatible with double.
   Sr_tmpv = Sr_tmp(mask);
   SXmat = repmat(Sr_tmpv.*xv,[1,Nt]);
   SYmat = repmat(Sr_tmpv.*yv,[1,Nt]);
   SZmat = repmat(Sr_tmpv.*zv,[1,Nt]);
   Dbmat = repmat(b((ncoils-1)*Nt+1:(ncoils)*Nt).',[Np,1]);
   Jx = Jx + SXmat.*A.*Dbmat/nx;
   Jy = Jy + SYmat.*A.*Dbmat/ny;
   Jz = Jz + SZmat.*A.*Dbmat/nz;
end
J = [Jx, Jy, Jz];
% dkx = 2*real(1i*Db'*A'*X*Sr*W*e);
% dky = 2*real(1i*Db'*A'*Y*Sr*W*e);
% dkz = 2*real(1i*Db'*A'*Z*Sr*W*e);
%
% db = A'*(A*b-d(mask));

