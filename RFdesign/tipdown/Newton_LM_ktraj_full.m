%% Newton method full update
% full update to k, alternate update b
% Levenberg-Marquardt method; add lambda*I to J'J and adjust lambda 
function [k,cost] = Newton_LM_ktraj_full(k, d, sens,mask,weightIm,beta,niter, xyzRange)

fieldmap = zeros(size(d)); % fieldmap is not used because the visiting time is known. 
tt = zeros(size(k,1),1)'; 
W = diag_sp(weightIm(mask)); 

b = 0; 
alp = 0.0001; 
alp0 = 1; 
% [dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap, tt, xyzRange);
[e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
dk = 2*real(1i*J'*W*e); 

grad_norm0 = norm([dk]); 
grad_norm = grad_norm0; 
b = dpls(As,d,beta,mask); 
e0 = norm(e); 
if ~isvar('niter')
    niter = 20; 
end


printf('iter | e | alp | norm(g)')
dobreak = 0; 
lam = 0.01; 
for iter = 1:niter
   alp = alp0; 
   % get the Hessian matrix
%    [H] = getHessian(k,b,d,A,e, sens,mask,weightIm);
%    [dkx, dky, dkz, db, e, A, J] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt,xyzRange);

   [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
   dk = 2*real(1i*J'*W*e); 
   H = real(J'*J);
   
   b = dpls(As,d,beta,mask);
   k_old = k; 
   kv_old = k_old(:); 
   e_old = norm(e); 
   grad_norm_old = norm([dk]); 

   while ((norm(e) == e_old) || (norm(e) > e_old)) 
      % || (grad_norm > grad_norm_old || grad_norm == grad_norm_old) 
   
   dire = (H + lam*eye(size(H)))\dk; %descending direction
   
   % update k
   kv = kv_old - alp * dire;
   k = reshape(kv,size(k)); 

   %    b = b -alp * db; 
   
   % evaluate the function
%    [dkxEva, dkyEva, dkzEva, dbEva, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
%    [dkxEva, dkyEva, dkzEva, dbEva, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap, tt, xyzRange);
   [e, A, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
   dkEva = 2*real(1i*J'*W*e); 

   grad_norm = norm([dkEva]); 
   if (norm(e) == e_old) || (norm(e) > e_old)
       lam = lam*10; % adjust lambda if e is larger or equal 
       lam
   end
   
   maxgo = max(alp * dire); 
   if maxgo < 0.00001
      dobreak = 1; 
      break; 
   end
   end % end of while 
   
   cost(iter) = norm(e); 
   printf('%d | %.5g | %.5g | %.3g \n', iter, cost(iter), lam, grad_norm) 
   lam = lam/10; 

   if grad_norm/grad_norm0 < 0.05
      break; 
   end
   if dobreak == 1; 
      break; 
   end
end 

cost = [e0 cost]; 


function [H] = getHessian(k,b,d,A,e,sens,mask,weightIm)
%% construct matrices
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

Db = diag_sp(b);
X = diag_sp(xx(mask));
Y = diag_sp(yy(mask)); 
Z = diag_sp(zz(mask)); 
Sr_tmp = sens(:,:,:,1); % working for single coil for now
Sr = diag_sp(Sr_tmp(mask)');
W = diag_sp(weightIm(mask));

%% construct Hessians
Hmiddle = A'*(X*Sr'*W*Sr*X)*A; 
% H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(conj(b).*diag(Hmiddle).*b)-2*real(full(e'*W*Sr*X*X*A*Db)');
% Hx = triu(H1,1)+tril(H1,-1)+diag(H2); 
H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(full(e'*W*Sr*X*X*A*Db)');
Hx = H1; % - diag(H2); 

Hmiddle = A'*(Y*Sr'*W*Sr*Y)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(full(e'*W*Sr*Y*Y*A*Db)');
Hy = H1; % - diag(H2); 

Hmiddle = A'*(Z*Sr'*W*Sr*Z)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(full(e'*W*Sr*Z*Z*A*Db)');
Hz = H1; % - diag(H2); 

Hmiddle = A'*(X*Sr'*W*Sr*Y)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(full(e'*W*Sr*X*Y*A*Db)');
Hxy = H1; % - diag(H2); 

Hmiddle = A'*(X*Sr'*W*Sr*Z)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(full(e'*W*Sr*X*Z*A*Db)');
Hxz = H1; % - diag(H2); 

Hmiddle = A'*(Z*Sr'*W*Sr*Y)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(full(e'*W*Sr*Z*Y*A*Db)');
Hzy = H1; % - diag(H2); 

H = [Hx, Hxy, Hxz; Hxy, Hy, Hzy; Hxz, Hzy, Hz]; 