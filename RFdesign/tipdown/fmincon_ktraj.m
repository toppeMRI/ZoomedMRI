%% Update k with constrains using fmincon
% full update to k (i.e., all kx, ky, kz), alternate update b
function [k,cost,nrmse] = fmincon_ktraj(k, d, sens,mask,weightIm,beta, niter, fieldmap, ndt, tt, xfov, yfov, zfov)
% ndt: interval between near k points
b = 0;
alp = 0.0001;
alp0 = 1;
[dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt);
grad_norm0 = norm([dkx; dky; dkz; db]);
grad_norm = grad_norm0;
b = dpls(A,d,beta,mask);
e0 = norm(e);

% the structure of finite diff C and bound from large to small are:
% first x,y,z, then positive/negative, then 1st/2nd order
nk = length(b);
C1 = Cdiff1(nk,'order',1);
C2 = Cdiff1(nk,'order',2);
Call = double(full([C1; C2]));
Call = [Call; -Call];
Call = [Call,zeros(size(Call)),zeros(size(Call)); zeros(size(Call)), Call, zeros(size(Call)); zeros(size(Call)), zeros(size(Call)), Call];
% Call = sparse(Call);

gamma = 4257.6; % Hz/Gauss
Gmax = 3.9; % Gauss/cm
Smax = 14.9e3; % Gauss/cm/s
dt = 4e-6; %s
bound1 = Gmax*gamma*dt*ndt; %*[ones(nk,1)*xfov,ones(nk,1)*yfov,ones(nk,1)*zfov]; % cycle/fov
bound2 = Smax*gamma*(dt)^2*ndt^2; %*[ones(nk,1)*xfov,ones(nk,1)*yfov,ones(nk,1)*zfov];
boundX = [bound1*ones(nk,1)*xfov; bound2*ones(nk,1)*xfov; bound1*ones(nk,1)*xfov; bound2*ones(nk,1)*xfov];
boundY = [bound1*ones(nk,1)*yfov; bound2*ones(nk,1)*yfov; bound1*ones(nk,1)*yfov; bound2*ones(nk,1)*yfov];
boundZ = [bound1*ones(nk,1)*zfov; bound2*ones(nk,1)*zfov; bound1*ones(nk,1)*zfov; bound2*ones(nk,1)*zfov];
bound = [boundX; boundY; boundZ];

% bound = sparse(bound(:));
% forwc = @(arg, x) x(end);
% backc = @(arg, y) zeros(size;
% C = fatrix2('forw', forwc, 'back', backc, 'idim', N, 'odim', N);

GetLast = zeros(1,3*nk);
GetLast = repmat(GetLast, [3,1]);
GetLast(1,nk) = 1;
GetLast(2,2*nk) = 1;
GetLast(3,3*nk) = 1;
% GetLast = sparse(GetLast);

kend = [0;0;0];
options = optimset('Algorithm','interior-point','Display','iter','MaxIter',50,...
   'GradObj','on', 'Hessian','user-supplied','HessFcn',@getHessian); % run interior-point algorithm
options = optimset('Algorithm','interior-point','Display','iter','MaxIter',3,...
   'GradObj','on'); % run interior-point algorithm

% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%         'Display','iter', 'GradObj','on', 'Hessian','user-supplied','HessFcn',@getHessian); % run interior-point algorithm

% project coeff to feasible region first
quadOptions = optimset('Algorithm', 'active-set'); 
[kv,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(k(:))), -k(:), Call, 0.99*bound,[],[],[],[],[],quadOptions); 
k = reshape(kv, [length(b),3]);

printf('iter | e')
for iter = 1: niter
   iter
   betaV = beta; % somehere beta has the same name as a matlab build in function and is shadowed in getHessian.
   save hessianParams b d A sens mask weightIm betaV fieldmap tt
   kv = k(:); 
   [kv,cost(iter)] = fmincon(@(kv) myfun(kv,b,d,sens,mask,weightIm,betaV,fieldmap,tt), kv, Call,bound, GetLast,kend,[],[],[],options);
   %    k = fminunc(@(k) myfun(k,b,d,sens,mask,weightIm,beta), k, options);
   k = reshape(kv, [length(b),3]);

   [dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt);
   b = dpls(A,d,beta,mask);
   nrmse(iter) = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
   printf('%d | %.5g | %.5g | %.3g \n', iter, cost(iter),nrmse)
end

end

function [f, dk] = myfun(kv,b,d, sens,mask,weightIm,betaV, fieldmap,tt)
k = reshape(kv,[length(b),3]);
[dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,betaV,fieldmap,tt);
f = norm(e)^2 + betaV*norm(b)^2;
dk = [dkx; dky; dkz];
end

function H = getHessian(k,lambda)
%% construct matrices
load hessianParams;
k = reshape(k, [length(k)/3,3]);
% notice beta here is fixed...
[dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,betaV,fieldmap,tt);

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
end
