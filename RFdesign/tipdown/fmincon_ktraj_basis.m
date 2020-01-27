%% Update k with constrains using basis function and fmincon
% full update to k (i.e., all kx, ky, kz), alternate update b
function [k,cost,nrmse,exetime,k_allIter, b_allIter]  = fmincon_ktraj_basis(k, d, sens,mask,weightIm,beta, niter, fieldmap, ndt, tt, xfov, yfov, zfov, xyzRange)

%% down sample in the space domain
scaleD = 0.5; % parameter to down sample the space dimension
d = imresize3(d, scaleD);
sens = imresize3(sens, scaleD);
mask = imresize3(mask, scaleD);
weightIm = imresize3(weightIm, scaleD);
fieldmap = imresize3(fieldmap, scaleD);
xyzRange.x = xyzRange.x(1:round(1/scaleD):end);
xyzRange.y = xyzRange.y(1:round(1/scaleD):end);

%% construct basis functions:
nCycle = 100;
coefs = eye(nCycle);
iEnd = nCycle + 1;
ti = linspace(0,iEnd,size(k,1))';
tis = repmat(ti, [1,nCycle]);
basis = bspline_1d_synth(coefs,tis,'ending','zero','order',2); figure(2),plot(ti, basis);
% basis = double(basis);
B = kron(eye(3,3),basis);
B = sparse(B); 

coeff = basis\k;
figure,plot(k),hold on,plot(basis*coeff,'--');
k0 = k; 
coeff = coeff(:);
coeff0 = coeff; % just record this value



%% constraint
% the structure of finite diff C and bound from large to small are:
% first x,y,z, then positive/negative, then 1st/2nd order
nk = size(k,1);
C1 = Cdiff1(nk,'order',1);
C2 = Cdiff1(nk,'order',2);
C1_basis = C1 * basis;
C2_basis = C2 * basis;

[Y, maxIdx] = max(C1_basis);
[Y, minIdx] = min(C1_basis);
peakIdxC1 = union(maxIdx,minIdx);
peakIdxC2 = peakIdxC1 + round((peakIdxC1(2) - peakIdxC1(1))/3);
nps = length(peakIdxC1); % number of peaks
C1_basis_peak = C1_basis(peakIdxC1,:);
C2_basis_peak = C2_basis(peakIdxC2,:);
Call = [C1_basis_peak; C2_basis_peak];
Call = double([Call; -Call]);
Call = [Call,zeros(size(Call)),zeros(size(Call)); zeros(size(Call)), Call, zeros(size(Call)); zeros(size(Call)), zeros(size(Call)), Call];

gamma = 4257.6; % Hz/Gauss
Gmax = 3.9; % Gauss/cm
Smax = 14.9e3; % Gauss/cm/s
dt = 4e-6; %s
bound1 = Gmax*gamma*dt*ndt; %*[ones(nps,1)*xfov,ones(nps,1)*yfov,ones(nps,1)*zfov]; % cycle/fov
bound2 = Smax*gamma*(dt)^2*ndt^2; %*[ones(nps,1)*xfov,ones(nps,1)*yfov,ones(nps,1)*zfov];
boundX = [bound1*ones(nps,1)*xfov; bound2*ones(nps,1)*xfov; bound1*ones(nps,1)*xfov; bound2*ones(nps,1)*xfov];
boundY = [bound1*ones(nps,1)*yfov; bound2*ones(nps,1)*yfov; bound1*ones(nps,1)*yfov; bound2*ones(nps,1)*yfov];
boundZ = [bound1*ones(nps,1)*zfov; bound2*ones(nps,1)*zfov; bound1*ones(nps,1)*zfov; bound2*ones(nps,1)*zfov];
bound = [boundX; boundY; boundZ]; % positive and negative

%% initial
b = 0;
% [dkx, dky, dkz, db, e, As, J] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt,xyzRange);
% [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
A = formA(k, sens, mask, fieldmap, tt, xyzRange);
[e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);

JB = J*B;
Np = length(weightIm(mask));
W = spdiag(weightIm(mask),0,Np,Np);
dkB = 2*real(1i*JB'*W*e);

b = dpls(As,d,beta,mask);
b0 = b; 
cost0 = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
nrmse0 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));

%% optimization
% kend = [0;0;0];
% options = optimset('Algorithm','interior-point','Display','iter','MaxIter',5,...
%    'GradObj','on', 'Hessian','user-supplied','HessFcn',@getHessian); % run interior-point algorithm
% options = optimset('Algorithm','interior-point','Display','iter','MaxIter',5,...
%    'GradObj','on'); % run interior-point algorithm
options = optimset('Algorithm','active-set','Display','iter','MaxIter',3,...
   'GradObj','on'); % run interior-point algorithm


% project coeff to feasible region first
%quadOptions = optimset('Algorithm', 'active-set');
quadOptions = optimset('Algorithm', 'interior-point-convex');
[coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -coeff, Call, 0.99*bound,[],[],[],[],[],quadOptions);
k = basis*reshape(coeff,[length(coeff)/3,3]);
A = formA(k, sens, mask, fieldmap, tt, xyzRange);
betaV = beta; % somehere beta has the same name as a matlab build in function and is shadowed in getHessian.
betavec = sqrt(betaV)*ones(size(b));

%%%%% the following few lines record the error right after B-spline fitting
[e, As, ~] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
b = dpls(As,d,beta,mask);
cost1 = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
nrmse1 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
k1 = k; 
b1 = b; 


tic
printf('iter | e')
for iter = 1: niter
   %    save hessianParams b d As sens mask weightIm betaV fieldmap tt
   [coeff,cost(iter)] = fmincon(@(coeff) myfun(coeff,b,d,sens,mask,weightIm,betaV, fieldmap,tt,B,xyzRange), coeff, Call,bound, [],[],[],[],[],options);
   %    k = fminunc(@(k) myfun(k,b,d,sens,mask,weightIm,beta), k, options);
   kv = B*coeff;
   k = reshape(kv, [length(b),3]);
   
   %    A = formA(k, sens, mask, fieldmap, tt, xyzRange);
   %    [e, As] = eval_func(A, b, d, sens, mask);
   
   [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
   JB = J*B;
   dkB = 2*real(1i*JB'*W*e);
   
   b = qpwls_pcg1_hao(b, As, W, d(mask), diag_sp(betavec),'niter',20);
   e = d(mask) - As*b; % record the residule error after update b.
   
   k_allIter{iter+2} = k;
   b_allIter{iter+2} = b; 
   nrmse(iter) = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
   costRec(iter) = norm(sqrt(W)*e)^2 + betaV*norm(b)^2; 
   printf('%d | %.5g | %.5g | %.3g \n', iter, costRec(iter), nrmse)
   
   exetime(iter) = toc;
end
k_allIter{1} = k0;
b_allIter{1} = b0;
k_allIter{2} = k1;
b_allIter{2} = b1;
costRec = [cost0 cost1 costRec]; 
nrmse = [nrmse0 nrmse1 nrmse]; 
exetime = [0 0 exetime]; 
save fmincon_exetime exetime nrmse cost

end

function [f, dkB] = myfun(coeff,b,d, sens,mask,weightIm,betaV, fieldmap,tt, B,xyzRange)
kv = B*coeff;
k = reshape(kv,[length(b),3]);
[e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask, weightIm, betaV, fieldmap, tt, xyzRange);
JB = J*B;
Np = length(weightIm(mask));
W = spdiag(weightIm(mask),0,Np,Np);
dkB = 2*real(1i*JB'*W*e);
f = norm(sqrt(W)*e)^2 + betaV*norm(b)^2;
end

function H = getHessian(k,lambda)
%% construct matrices
% load hessianParams;
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
