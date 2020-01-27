%% Newton method full update (precondition gradient descent with backtracking)
% full update to k, alternate update b
% Levenberg-Marquardt method; add lambda*I to hessian and adjust lambda

function [k,costA, nrmse] = Newton_LM_ktraj_basis(k, d, sens,mask,weightIm,beta,niter,fieldmap,tt, xyzRange)

% k = zeros(size(k));

% construct basis functions:
if 0 % shiftted gaussian
   kLen = size(k,1);
   nOneCycle = 20;
   centers = [0.5:(kLen/nOneCycle)-0.5]*nOneCycle;
   nbasis = length(centers);
   ki = 1:kLen;
   width = 60;
   kiMat = repmat(ki',[1,nbasis]);
   centersMat = repmat(centers,[kLen,1]);
   basis = exp(-(kiMat-centersMat).^2/width);
   
   coeff0 = basis\k;
   figure,plot(k),hold on,plot(basis*coeff0,'--');
   
   B = kron(eye(3,3),basis);
   k = basis*coeff0;
   
else % bspline
   nCycle = 50;
   coefs = eye(nCycle);
   iEnd = nCycle + 1;
   ti = linspace(0,iEnd,size(k,1))';
   tis = repmat(ti, [1,nCycle]);
   basis = bspline_1d_synth(coefs,tis,'ending','zero','order',2); figure(2),plot(ti, basis);
   basis = double(basis); 
   B = kron(eye(3,3),basis);
   
   coeff0 = basis\k;
   figure,plot(k),hold on,plot(basis*coeff0,'--');
   k = basis*coeff0;
end

b = 0;
[dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt,xyzRange);
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
beta2 = 0;
cost0 = norm(e)^2 + beta*norm(b)^2;
nrmse0 = norm(e)/sqrt(numel(e))/max(abs(d(mask))); 

printf('iter | cost | nrmse | lam | norm(g)')
dobreak = 0;
lam = 1e-10;
coeff = coeff0(:);
alp0 = 1;
alp = alp0;

for iter = 1:niter
   % get the Hessian matrix for k
   kv = k(:);
   [H] = getHessian(k,b,d,A,e, sens,mask,weightIm);
   HB = B'*H*B;
   
   [dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
   dkB = B'*[dkx; dky; dkz];
   
   coeff_old = coeff;
   
   cost = norm(e)^2 + beta*norm(b)^2;
   cost_old = cost;
   grad_norm_old = norm(dkB);
   
   while (cost == cost_old) || (cost > cost_old)
      % || (grad_norm > grad_norm_old || grad_norm == grad_norm_old)
      
      dire = (HB + lam*eye(size(HB)))\dkB; %descending direction
      
      % update coeffs
      coeff = coeff_old - alp * dire;
     
      kv = B*coeff;
      k = reshape(kv,size(k));
      
      % evaluate the function
      [dkxEva, dkyEva, dkzEva, dbEva, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
      grad_norm = norm([dkxEva; dkyEva; dkzEva]);
      cost = norm(e)^2 + beta*norm(b)^2;
      
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
   
   b = dpls(A,d,beta,mask); %update b
   e = d(mask) - A*b; % record the residule error after update b. 

   %%
   costA(iter) = cost;
   nrmse(iter) = norm(e)/sqrt(numel(e))/max(abs(d(mask))); 
   printf('%d | %.5g | %.5g | %.5g | %.3g \n', iter, costA(iter), nrmse(iter), lam, grad_norm)
   lam = lam/10;
   
   if grad_norm/grad_norm0 < 0.05
      break;
   end
   if dobreak == 1;
      break;
   end
end

costA = [cost0 costA];
nrmse = [nrmse0 nrmse]; 
save testA A b d k e beta mask


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
% H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(conj(b).*diag(Hmiddle).*b)-2*real(full(e'*W*Sr*X*X*A*Db)');
% Hx = triu(H1,1)+tril(H1,-1)+diag(H2);
% Hx = H1; % - diag(H2);
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
