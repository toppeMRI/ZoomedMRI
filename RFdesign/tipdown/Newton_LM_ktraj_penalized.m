%% Newton method full update (precondition gradient descent with backtracking)
% full update to k, alternate update b
% Levenberg-Marquardt method; add lambda*I to hessian and adjust lambda

function [k,costA] = Newton_LM_ktraj_penalized(k, d, sens,mask,weightIm,beta,niter)
b = 0; 
alp = 0.0001; 
alp0 = 1; 
[dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm,beta);
grad_norm0 = norm([dkx; dky; dkz; db]); 
grad_norm = grad_norm0; 
b = dpls(A,d,beta,mask); 
if ~isvar('niter')
    niter = 20; 
end


kLen = size(k,1); 
I = speye(kLen,kLen);
C1 = I - circshift(I, [1,0]);
C2 = C1*C1; 
D1 = kron(speye(3,3),C1); 
D2 = kron(speye(3,3),C2); 
beta1 = 2; 
beta2 = 1; 
cost0 = norm(e) + beta2*norm(D2*k(:)); 

printf('iter | e | alp | norm(g)')
dobreak = 0; 
lam = 0.01; 
for iter = 1:niter
   alp = alp0; 
   % get the Hessian matrix for k
   kv = k(:); 
   [H] = getHessian(k,b,d,A,e, sens,mask,weightIm);
   Ha = H + beta2*D2'*D2; 
   [dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
   dkxA = dkx+beta2*C2'*k(:,1); 
   dkyA = dky+beta2*C2'*k(:,2); 
   dkzA = dkz+beta2*C2'*k(:,3); 
   b = dpls(A,d,beta,mask); %update b
   kv_old = k(:); 
   
   cost = norm(e) + beta2*norm(D2*kv); 
   cost_old = cost; 
   grad_norm_old = norm([dkxA; dkyA; dkzA]); 

   while (cost == cost_old) || (cost > cost_old) 
      % || (grad_norm > grad_norm_old || grad_norm == grad_norm_old) 
   
   dire = (Ha+ lam*eye(size(Ha)))\[dkxA;dkyA;dkzA]; %descending direction
   
   % update k
   kv = kv_old - alp * dire; 
   k = reshape(kv,size(k)); 

   %    b = b -alp * db; 
   
   % evaluate the function
   [dkxEva, dkyEva, dkzEva, dbEva, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
   grad_norm = norm([dkxEva; dkyEva; dkzEva]); 
   cost = norm(e) + beta2*norm(D2*kv); 

   if (cost == cost_old) || (cost > cost_old) 
       lam = lam*10; % adjust lambda if e is larger or equal 
       lam
   end
   
   maxgo = max(alp * dire); 
   if maxgo < 0.00001
      dobreak = 1; 
      break; 
   end
   end % end of while 
   
 %%  
   costA(iter) = cost; 
   printf('%d | %.5g | %.5g | %.3g \n', iter, costA(iter), lam, grad_norm) 
   lam = lam/10; 

   if grad_norm/grad_norm0 < 0.05
      break; 
   end
   if dobreak == 1; 
      break; 
   end
end 

costA = [cost0 costA]; 


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
Sr_tmp = sens(:,:,:,1); 
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
