%% BBGM with nonmontone line search
% following paper by Marcos Raydan SIAM 1997
% alternatively update k and b
% not finished.

function [k,costA, nrmse] = BBGM_proj_ktraj_basis(k, d, sens,mask,weightIm,beta,niter,fieldmap,tt, xyzRange, ndt, xfov, yfov, zfov)

if ~isvar('niter')
   niter = 20;
end
% k = zeros(size(k));

% construct basis functions:
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

% construct difference matrix and bound
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
[dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt,xyzRange);
grad_norm0 = norm([dkx; dky; dkz; db]);
grad_norm = grad_norm0;
b = dpls(A,d,beta,mask);
cost0 = norm(e)^2 + beta*norm(b)^2;
nrmse0 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));


%% projected BBGM
printf('iter | cost | nrmse | lam | norm(g)')
dobreak = 0;
lam = 1e-10;
coeff = coeff0(:);
alp0 = 1;
alp = alp0;

quadOptions = optimset('Algorithm', 'active-set');
[coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -coeff, Call, 0.99*bound,[],[],[],[],[],quadOptions);
k = basis*reshape(coeff,[length(coeff)/3,3]);

scale = 1.2; % table is constructed from scale*kmin to scale*kmax
dimTable = 8;
[HessianMiddleTable, kspacing] = createHessianTable(k,b,sens,mask,weightIm,fieldmap, tt, xyzRange, scale, dimTable);

% [H] = getHessian(k,b,d,A,e, sens,mask,weightIm, xyzRange);
% HB = B'*H*B;



for iter = 1:niter
   % get the Hessian matrix for k
   kv = k(:);
   %     [H] = getHessian(k,b,d,A,e, sens,mask,weightIm, xyzRange);
   %     [H] = lookUpHessian(k,b,HessianMiddleTable,kspacing,scale,dimTable);
   %     HB = B'*H*B;
   
   % calculate gradient
   [dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
   dkB = B'*[dkx; dky; dkz];
   g = dkB;
   coeff_old = coeff;
   
   cost = norm(e)^2 + beta*norm(b)^2;
   grad_norm = norm(g);
   
   
   % update of k using BBGM

   % BBGM parameters
   kgrad_thre = 5;
   alpha = 10; gam = 1e-4; eps = 1e-10; M = 10; sig1 = 0.1; sig2 = 0.5;
   iter_bbgm = 0;
   cost_all = [] ; % for bbgm iterations
   while (grad_norm > kgrad_thre)
      iter_bbgm = iter_bbgm+1;
      cost_all(iter_bbgm) = cost; 
      
      % reset alpha if needed
      if (alpha < eps || alpha > 1/eps)
         alpha = eps;
      end
      
      % step size
      lam = 1/alpha;
      
      % nonmonotone line search
      while(true)
         
         % update for candidate
         coeff = coeff_old - lam * g;
         
         % projection
%          [coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -coeff, Call, 0.99*bound,[],[],[],[],[],quadOptions);
         kv = B*coeff;
         k = reshape(kv, size(k));
         
         % evaluate the function
         [dkx, dky, dkz, db, e, A] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);      
         cost = norm(e)^2 + beta*norm(b)^2;
         grad_norm = norm(g);
         
         % check the search criteria
         if(cost < (max(cost_all(max([1,iter_bbgm - M]):end)) - gam*lam*(g'*g)))
            g_old = g; % update g 
            dkB = B'*[dkx; dky; dkz];
            g = dkB;
            coeff_old = coeff; 
            alpha = -(g_old'*(g - g_old) / (lam*(g_old'*g_old)));
            break; 
         else
            sig = mean([sig1 sig2]);
            lam = sig*lam;
         end
         
      end
      printf('%d | %.5g | %.5g | %.5g | \n', iter_bbgm, cost_all(iter_bbgm), lam, grad_norm)

   end
   
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
