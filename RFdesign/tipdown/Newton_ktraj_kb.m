%% Newton method (precondition gradient descent with backtracking)
function [k,cost] = PGD_ktraj_kb(k, d, sens, mask, weightIm)
b = 0; 
alp = 0.0001; 
alp0 = 2; 
[dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm); 
grad_norm0 = norm([dkx; dky; dkz; db]); 
grad_norm = grad_norm0; 
b = pinv(A)*d(mask); 
e0 = norm(e); 

printf('iter | e | alp | norm(g)')
dobreak = 0; 
for iter = 1:50
   alp = alp0; 
   % get the Hessian matrix
   [Hx, Hy, Hz, Hb] = getHessian(k,b,d,A,e, sens,mask,weightIm); 
   [dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm);

%    b = pinv(A)*d(mask);
   k_old = k; 
   b_old = b; 
   e_old = norm(e); 
   grad_norm_old = grad_norm; 
   while ((norm(e) == e_old) || (norm(e) > e_old))
         % || (grad_norm > grad_norm_old || grad_norm == grad_norm_old) 

   % backtracking step size
   
   k(:,1) = k_old(:,1) - alp * (Hx\dkx);
   k(:,2) = k_old(:,2) - alp * (Hy\dky); 
   k(:,3) = k_old(:,3) - alp * (Hz\dkz); 
   b = b_old - alp * (Hb\db); 
   % get the gradient 
   [dkxEva, dkyEva, dkzEva, dbEva, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
   grad_norm = norm([dkxEva; dkyEva; dkzEva; dbEva]); 
   alp = alp/2;
   maxgo = max([alp * (Hx\dkx); alp*(Hy\dky); alp*(Hz\dkz)]); 
   if maxgo < 0.001
      dobreak = 1; 
      printf('too small step'); 
      break; 
   end
   end
   
   cost(iter) = norm(e); 
   printf('%d | %.5g | %.5g | %.3g \n', iter, cost(iter), alp*2, grad_norm) 

   if grad_norm/grad_norm0 < 0.05
      break; 
      printf('small enough gradient'); 
   end
   if dobreak == 1; 
      break; 
   end
end 

cost = [e0 cost]; 

function [Hx, Hy, Hz,Hb] = getHessian(k,b,d,A,e,sens,mask,weightIm)
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
H2 = 2*real(full(e'*W*Sr*X*X*A*Db)');
Hx = H1 - diag(H2); 

Hmiddle = A'*(Y*Sr'*W*Sr*Y)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
H2 = 2*real(full(e'*W*Sr*Y*Y*A*Db)');
Hy = H1 - diag(H2); 

Hmiddle = A'*(Z*Sr'*W*Sr*Z)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
H2 = 2*real(full(e'*W*Sr*Z*Z*A*Db)');
Hz = H1 - diag(H2); 

Hb = A'*A; 

