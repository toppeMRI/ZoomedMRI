function [e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta,fieldmap,tt, xyzRange)
% calculate gradient over kx, ky, kz, and b
% function [dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt, xyzRange)

%start=toc;

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


ncoils = size(sens,4);
Nt = size(A,2); 
%% construct the matrix A and evaluate residual
% combine sens and A into As
sensi = sens(:,:,:,1);
Atmp = spdiag(sensi(mask))*A;
for ic = 2:ncoils
   sensi = sens(:,:,:,ic);
   Atmp = [Atmp; spdiag(sensi(mask))*A];
end
As = reshape(Atmp, [sum(mask(:)),ncoils*Nt]);

if b == 0  % calculate an initialization for b if set b = 0 as the input .
   b = double(dpls(As,single(d),beta,mask));
end

% calculate error
% a = embed(A*b, mask);
% r = d - embed(A*b, mask);
% e = r(mask);
e = d(mask) - As*b;

%% Compute the gradients and Jacobian
Np = length(xx(mask));

xv = xx(mask)/nx;
yv = yy(mask)/ny;
zv = zz(mask)/nz;

%X = spdiag(xx(mask)/nx,0,Np,Np);
%Y = spdiag(yy(mask)/ny,0,Np,Np);
%Z = spdiag(zz(mask)/nz,0,Np,Np);
% W = spdiag(weightIm(mask),0,Np,Np);

% These are the 'SXAB' in Eq. (10), i.e., it doesn't include the last multiplication by H
Jx = zeros(Np, Nt);   % Np = # spatial locations in mask;  Nt = number of time-samples in b
Jy = Jx;
Jz = Jx;

%fprintf(1,'inside kb_Jacobian_inputA(): %f sec\n', toc-start);

for ic = 1:ncoils
   %Db = spdiag(b((ncoils-1)*Nt+1:(ncoils)*Nt),0,Nt,Nt);
   %Db = diag_sp(b((ncoils-1)*Nt+1:(ncoils)*Nt));   % Also works
   Sr_tmp = sens(:,:,:,ic);
   
   % % sparse matrix only works with double
%       Sr = spdiag(Sr_tmp(mask),0,Np,Np);
%       Jx = Jx + Sr*X*double(A)*Db; % remove the negative sign here
%       Jy = Jy + Sr*Y*double(A)*Db;
%       Jz = Jz + Sr*Z*double(A)*Db;
   
   % % avoid sparse matrix because it is not compatible with double.
	% the following is faster than the above
   Sr_tmpv = Sr_tmp(mask);
   SXmat = repmat(Sr_tmpv.*xv,[1,Nt]);
   SYmat = repmat(Sr_tmpv.*yv,[1,Nt]);
   SZmat = repmat(Sr_tmpv.*zv,[1,Nt]);
   Dbmat = repmat(b((ncoils-1)*Nt+1:(ncoils)*Nt).',[Np,1]);
	%times_inplace_mex(A,Dbmat);   % Dbmat = A.*Dbmat
   %Jx = Jx + SXmat.*(A.*Dbmat);
   %Jy = Jy + SYmat.*(A.*Dbmat);
   %Jz = Jz + SZmat.*(A.*Dbmat);
	Dbmat = A.*Dbmat; % in-place calculation 
   Jx = Jx + SXmat.*(Dbmat);
   Jy = Jy + SYmat.*(Dbmat);
   Jz = Jz + SZmat.*(Dbmat);
end
%start = toc;
J = [Jx, Jy, Jz];
%fprintf(1,'inside kb_Jacobian_inputA(): %f sec\n', toc-start);
%J = zeros(Np,3*Nt);
%J(:,1:Nt) = Jx;
%J(:,(Nt+1):(2*Nt)) = Jy;
%J(:,(2*Nt+1):(3*Nt)) = Jz;
%fprintf(1,'inside kb_Jacobian_inputA(): %f sec\n', toc-start);

%fprintf(1,'inside kb_Jacobian_inputA(): %f sec\n', toc-start);
% dkx = 2*real(1i*Db'*A'*X*Sr*W*e);
% dky = 2*real(1i*Db'*A'*Y*Sr*W*e);
% dkz = 2*real(1i*Db'*A'*Z*Sr*W*e);
%
% db = A'*(A*b-d(mask));

