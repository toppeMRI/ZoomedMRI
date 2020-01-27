%% k update in the ADMM 
% modified from Newton_LM_ktraj_full. 
% Newton method full update (precondition gradient descent with backtracking)
% full update to k, alternate update b
% Levenberg-Marquardt method; add lambda*I to hessian and adjust lambda 
function [k,cost] = ADMM_updateK_ktraj(k, d, sens,mask,weightIm,beta,lambda,mu,C,xfov,yfov,zfov)

b = 0; 
alp = 0.0001; 
alp0 = 1; 
[dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm,beta);
grad_norm0 = norm([dkx; dky; dkz; db]); 
grad_norm = grad_norm0; 
b = dpls(A,d,beta,mask); 
e0 = norm(e); 

Glimit = 4; %G/cm
Slimit = 15; %G/cm/ms
dt = 4*4e-3; %ms
gamma = 4257; %Hz/G 
kdmax = Glimit*gamma*dt/1000; % cycle/cm 
kddmax = Slimit*gamma*dt*dt/1000; %cycle/cm
kdmaxV = [kdmax*ones(length(k),1)*xfov; kdmax*ones(length(k),1)*yfov; kdmax*ones(length(k),1)*zfov];
kddmaxV = [kddmax*ones(length(k),1)*xfov; kddmax*ones(length(k),1)*yfov; kddmax*ones(length(k),1)*zfov];
ind = C*k(:) - lambda/mu > 0; % if (Ck)_i - lambda/mu > 0; 
v = C*k(:) - [kdmaxV; kddmaxV]-lambda/mu; 
v(ind) = 0; 
dkp = mu*C'*v;
hkp = mu*full(C'*C); 
printf('iter | e | alp | norm(g)')
dobreak = 0; 
lam = 0.01; 
for iter = 1:5
   alp = alp0; 
   % get the Hessian matrix
   [H] = getHessian(k,b,d,A,e, sens,mask,weightIm);
   H = H  + full(hkp); 
   [dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
   b = dpls(A,d,beta,mask);
   k_old = k; 
   kv_old = k_old(:); 
   e_old = norm(e); 
   grad_norm_old = norm([dkx; dky; dkz; db]); 

   while ((norm(e) == e_old) || (norm(e) > e_old)) 
      % || (grad_norm > grad_norm_old || grad_norm == grad_norm_old) 
   
   dire = (H+ lam*eye(size(H)))\([dkx;dky;dkz]+dkp); %descending direction
   
   % update k
   kv = kv_old - alp * dire;
   k = reshape(kv,size(k)); 

   %    b = b -alp * db; 
   
   % evaluate the function
   [dkxEva, dkyEva, dkzEva, dbEva, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
   grad_norm = norm([dkxEva; dkyEva; dkzEva; dbEva]); 
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