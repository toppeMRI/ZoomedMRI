% compare the analytical gradient to the numerical gradient
% function [HdiagNum, Hdiag] = test_hessian(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
nInd = 150;
for i = 1:nInd
   ind = i;
   dd = 1e-8;
   coeff1 = coeff; 
   coeff2 = coeff1; 
   coeff2(ind) = coeff2(ind) + dd; 
   
   kv = B*coeff1;
   k1 = reshape(kv,size(k));
   
   kv = B*coeff2; 
   k2 = reshape(kv,size(k));
   
   W = diag_sp(weightIm(mask));
   
   [e1, As, J1] = kb_Jacobian_withB0(k1, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
   JB1 = J1*B;
   dkB1 = 2*(2*pi)*real(1i*JB1'*W*e1);

   [e2, As, J2] = kb_Jacobian_withB0(k2, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
   JB2 = J2*B;
   dkB2 = 2*(2*pi)*real(1i*JB2'*W*e2);

   HdiagNum(i) = (dkB2(i) - dkB1(i))/dd; 
end
   HB = 2*(2*pi)^2*real(JB1'*W*JB1);   
   Hdiag = diag(HB); 
   
figure,plot(HdiagNum); hold on; plot(Hdiag,'r')
% figure,plot(Hdiag./HdiagNum); 


% function [H] = getHessian(k,b,d,A,e,sens,mask,weightIm,xyzRange)
%% construct matrices
nx = size(d,1);
ny = size(d,2);
nz = size(d,3);
x = xyzRange.x/(xyzRange.x(2)-xyzRange.x(1))/nx;
y = xyzRange.y/(xyzRange.y(2)-xyzRange.y(1))/ny;
z = xyzRange.z/(xyzRange.z(2)-xyzRange.z(1))/nz;
[xx,yy,zz]=ndgrid(x,y,z);

kx = k(:,1);
ky = k(:,2);
kz = k(:,3);

Db = spdiag(b);
X = spdiag(xx(mask));
Y = spdiag(yy(mask));
Z = spdiag(zz(mask));
Sr_tmp = sens(:,:,:,1);
Sr = spdiag(Sr_tmp(mask)');
W = spdiag(weightIm(mask));

%% construct Hessians
SWS = Sr'*W*Sr;
Hmiddle = A'*(X*SWS*X)*A;
H1 = 2*(2*pi)^2*real(full(Db'*Hmiddle*Db));
H2 = 2*(2*pi)^2*real(conj(b).*diag(Hmiddle).*b)-2*(2*pi)^2*real(full(e'*W*Sr*X*X*A*Db)');
Hx = triu(H1,1)+tril(H1,-1)+diag(H2);
% Hx = H1 - diag(H2);
HxB = B(1:end/3,1:end/3)'*Hx*B(1:end/3,1:end/3); 
Hx = 2*real(Db'*Hmiddle*Db);

Hmiddle = A'*(Y*SWS*Y)*A;
Hy = 2*real(Db'*Hmiddle*Db);

Hmiddle = A'*(Z*SWS*Z)*A;
Hz = 2*real(Db'*Hmiddle*Db);

Hmiddle = A'*(X*SWS*Y)*A;
Hxy = 2*real(Db'*Hmiddle*Db);

Hmiddle = A'*(X*SWS*Z)*A;
Hxz = 2*real(Db'*Hmiddle*Db);

Hmiddle = A'*(Z*SWS*Y)*A;
Hzy = 2*real(Db'*Hmiddle*Db);

H = [Hx, Hxy, Hxz; Hxy, Hy, Hzy; Hxz, Hzy, Hz];
