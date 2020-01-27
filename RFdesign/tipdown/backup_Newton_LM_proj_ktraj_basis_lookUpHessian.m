%% Projected LM 
% express k-traj using b-spline basis
% Levenberg-Marquardt method; add lambda*I to J'J and adjust lambda

function [k,costA, nrmse] = backup_Newton_LM_proj_ktraj_basis_lookUpHessian(k, d, sens,mask,weightIm,beta,niter,fieldmap,tt, xyzRange, ndt, xfov, yfov, zfov)

if ~isvar('niter')
   niter = 20;
end
% k = zeros(size(k));



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
nCycle = 50;
coefs = eye(nCycle);
iEnd = nCycle + 1;
ti = linspace(0,iEnd,size(k,1))';
tis = repmat(ti, [1,nCycle]);
basis = bspline_1d_synth(coefs,tis,'ending','zero','order',2); figure(2),plot(ti, basis);
% basis = single(basis); 
B = kron(eye(3,3),basis);

coeff0 = basis\k;
figure,plot(k),hold on,plot(basis*coeff0,'--');
k = basis*coeff0;

%% construct difference matrix and bound
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


%% get initial
b = 0;
% [dkx, dky, dkz, db, e, As, J] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt,xyzRange);
% [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
A = formA(k, sens, mask, fieldmap, tt, xyzRange); 
[e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);       

JB = J*B;
W = diag_sp(weightIm(mask)); 
dkB = 2*real(1i*JB'*W*e); 
      
grad_norm0 = norm(dkB);
grad_norm = grad_norm0;
b = dpls(As,d,beta,mask);
cost0 = norm(e)^2 + beta*norm(b)^2;
nrmse0 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));

betavec = sqrt(beta)*ones(size(b));
%% projected LM
printf('iter | cost | nrmse | lam | norm(g)')
dobreak = 0;
lam = 1e-10;
coeff = coeff0(:);
alp0 = 1;
alp = alp0;

quadOptions = optimset('Algorithm', 'interior-point-convex','Display','off');
[coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -double(coeff), Call, double(0.99*bound),[],[],[],[],[],quadOptions);
k = basis*reshape(coeff,[length(coeff)/3,3]);
 
exetime = 0;
tic; 
niterK = 5;
for iter = 1:niter
   
   % optimize k
   
%    [H] = getHessian(k,b,d,A,e, sens,mask,weightIm, xyzRange);
%    [H] = lookUpHessian(k,b,HessianMiddleTable,kspacing,scale,dimTable);
%    HB = B'*H*B;
   
%    momentum = 2;
%    coeff_z = coeff; 
   for iterK = 1:niterK
      % get the Hessian matrix for k

%       [dkx, dky, dkz, db, e, A, J] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
%       dkB = B'*[dkx; dky; dkz];
         
%       [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);    
      [e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);       
      JB = J*B; 
      dkB = 2*(2*pi)*real(1i*JB'*W*e); 
      HB = 2*(2*pi)^2*real(JB'*W*JB);
      
%       L = eigs(HB,1)*speye(size(HB)); 
      
      coeff_old = coeff;
      
      cost = norm(e)^2 + beta*norm(b)^2;
      cost_old = cost;
      grad_norm_old = norm(dkB);
      
%       momentum_old = momentum; % nesterov momentum
%       momentum = (1 + sqrt(1 + 4*momentum_old^2))/2; 
%       coeff_z_old = coeff_z; 
      
      while (cost == cost_old) || (cost > cost_old)
         % || (grad_norm > grad_norm_old || grad_norm == grad_norm_old)
         
         dire = (HB + lam*eye(size(HB)))\dkB; %descending direction
%          dire = (HB + lam*diag(diag(HB)))\dkB; %descending direction
%          dire = (L + lam*speye(size(HB)))\dkB; %descending direction
         
         % update coeffs
         coeff = coeff_old - alp * dire;
%          coeff = coeff_z_old - dire; 
%          coeff_z = coeff + (momentum_old - 1)/momentum*(coeff - coeff_old); 
         
         % projection
         [coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -double(coeff), Call, double(0.99*bound),[],[],[],[],[],quadOptions);
         
         
         kv = B*coeff;
         k = reshape(kv,size(k));
         
         % evaluate the function
         A = formA(k, sens, mask, fieldmap, tt, xyzRange); 
         [e, As] = eval_func(A, b, d, sens, mask);
         grad_norm = norm(dkB);
         cost = e'*W*e + beta*norm(b)^2;
         
         if (cost == cost_old) || (cost > cost_old)
            lam = lam*10; % adjust lambda if e is larger or equal
            lam
         end
         
         maxgo = max(abs(alp * dire));
         if maxgo < 0.00001
            dobreak = 1;
            warning('step size is too small'); 
            break;
         end
      end % end of while
      
   end % end of niterK operations

   %    b2 = dpls(A,d,beta,mask); %update b
   b = qpwls_pcg1_hao(b, As, W, d(mask), diag_sp(betavec),'niter',20); 
   e = d(mask) - As*b; % record the residule error after update b.
   
   %%
   costA(iter) = cost;
   nrmse(iter) = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
   exetime(iter) = toc; 
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
exetime = [0 exetime]; 
save projectd_LM_exetime costA nrmse exetime


function [H] = getHessian(k,b,d,A,e,sens,mask,weightIm,xyzRange)
%% construct matrices
nx = size(d,1);
ny = size(d,2);
nz = size(d,3);
x = xyzRange.x/(xyzRange.x(2)-xyzRange.x(1));
y = xyzRange.y/(xyzRange.y(2)-xyzRange.y(1));
z = xyzRange.z/(xyzRange.z(2)-xyzRange.z(1));
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

function [HessianMiddleTable, kspacing] = createHessianTable(k,b,sens,mask,weightIm,fieldmap, tt, xyzRange, scale, dimTable)
%% construct matrices
% nx = size(d,1);
% ny = size(d,2);
% nz = size(d,3);
% x = [ceil(-nx/2):1:floor((nx-1)/2)];
% y = [ceil(-ny/2):1:floor((ny-1)/2)];
% z = [ceil(-nz/2):1:floor((nz-1)/2)];

x = xyzRange.x/(xyzRange.x(2)-xyzRange.x(1));
y = xyzRange.y/(xyzRange.y(2)-xyzRange.y(1));
z = xyzRange.z/(xyzRange.z(2)-xyzRange.z(1));
[xx,yy,zz]=ndgrid(x,y,z);

kx = k(:,1);
ky = k(:,2);
kz = k(:,3);


kxMax = scale*max(kx);
kxMin = scale*min(kx);
kyMax = scale*max(ky);
kyMin = scale*min(ky);
kzMax = scale*max(kz);
kzMin = scale*min(kz);


kxTable = linspace(kxMin, kxMax, dimTable);
kyTable = linspace(kyMin, kyMax, dimTable);
kzTable = linspace(kzMin, kzMax, dimTable);
kspacing = [kxTable(2) - kxTable(1), kyTable(2) - kyTable(1), kzTable(2) - kzTable(1)];

[kxx, kyy, kzz] = ndgrid(kxTable, kyTable, kzTable);
kk = [kxx(:), kyy(:), kzz(:)];
A = formA(kk, sens, mask, zeros(size(fieldmap)), zeros(size(kxx(:)))', xyzRange);

Db = spdiag(b);
X = spdiag(xx(mask));
Y = spdiag(yy(mask));
Z = spdiag(zz(mask));
Sr_tmp = sens(:,:,:,1);
Sr = spdiag(Sr_tmp(mask)');
W = spdiag(weightIm(mask));

%% construct Hessians table
SWS = Sr'*W*Sr;
Hxm = A'*(X*SWS*X)*A;
% H1 = 2*real(full(Db'*Hmiddle*Db));
% H2 = 2*real(conj(b).*diag(Hmiddle).*b)-2*real(full(e'*W*Sr*X*X*A*Db)');
% Hx = triu(H1,1)+tril(H1,-1)+diag(H2);
% Hx = H1; % - diag(H2);

% Hx = 2*real(Db'*Hmiddle*Db);

Hym = A'*(Y*SWS*Y)*A;

Hzm = A'*(Z*SWS*Z)*A;

Hxym = A'*(X*SWS*Y)*A;

Hxzm = A'*(X*SWS*Z)*A;

Hzym = A'*(Z*SWS*Y)*A;

HessianMiddleTable = [Hxm, Hxym, Hxzm; Hxym, Hym, Hzym; Hxzm, Hzym, Hzm];

function [H] = lookUpHessian(k,b, HessianMiddleTable,kspacing,scale,dimTable)
kx = k(:,1);
ky = k(:,2);
kz = k(:,3);

kxMax = scale*max(kx);
kxMin = scale*min(kx);
kyMax = scale*max(ky);
kyMin = scale*min(ky);
kzMax = scale*max(kz);
kzMin = scale*min(kz);

kxIdx = round((kx - kxMin)/kspacing(1))+1;
kyIdx = round((ky - kyMin)/kspacing(2))+1;
kzIdx = round((kz - kzMin)/kspacing(3))+1;
sizeKgrid = size(HessianMiddleTable,1)/3;
kIdxOne = kxIdx + (kyIdx-1)*dimTable + (kzIdx-1)*dimTable^2;
kIdx = [kIdxOne; kIdxOne+sizeKgrid; kIdxOne+2*sizeKgrid];
Hmiddle = HessianMiddleTable(kIdx,kIdx);

Db = spdiag([b;b;b]);

H = 2*real(Db'*Hmiddle*Db);
