%% Solve a quadratic program at each iteration. Accelerate using momentum (FISTA). JFN 2015

function [k, costRec, nrmse, nrmse_ov, exetime, k_allIter, b_allIter] = fista_ktraj_basis(k, d, sens,mask,weightIm,beta,niter,fieldmap,tt, xyzRange, ndt, xfov, yfov, zfov)

tic   % needed when using parfor loop

load rois;    % iv_roi, ov_roi, edge_roi (created in main_ktCont.m)

bsplineorder = 2;

if ~isvar('niter')
   niter = 20;
end
% k = zeros(size(k));

k0 = k; 


%% down sample in the space domain
scaleD = 0.5; % parameter to down sample the space dimension
d = imresize3(d, scaleD);
sens = imresize3(sens, scaleD);
mask = imresize3(mask, scaleD);
weightIm = imresize3(weightIm, scaleD);
fieldmap = imresize3(fieldmap, scaleD);
xyzRange.x = xyzRange.x(1:round(1/scaleD):end);
xyzRange.y = xyzRange.y(1:round(1/scaleD):end);
iv_roi = imresize3(iv_roi,scaleD);
ov_roi = imresize3(ov_roi,scaleD);
edge_roi = imresize3(edge_roi,scaleD);



%% construct basis functions:
if bsplineorder==2
	nCycle = 50;
else
	nCycle = 20;
end
coefs = eye(nCycle);
iEnd = nCycle + 1;
ti = linspace(0,iEnd,size(k,1))';
tis = repmat(ti, [1,nCycle]);
basis = bspline_1d_synth(coefs,tis,'ending','zero','order',bsplineorder); %figure(2),plot(ti, basis);
% basis = double(basis);
B = kron(eye(3,3),basis);
B = sparse(double(B)); 

coeff0 = basis\k;
%figure,plot(k),hold on,plot(basis*coeff0,'--'); title('fitting to init-k'); 

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
if bsplineorder==2
	nps = length(peakIdxC1); % number of peaks
else
	nps = size(basis,1);
end
C1_basis_peak = C1_basis(peakIdxC1,:);
C2_basis_peak = C2_basis(peakIdxC2,:);
if bsplineorder==2
	Call = [C1_basis_peak; C2_basis_peak];
else
	Call = [C1_basis; C2_basis];
end
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
A = formA(k0, sens, mask, fieldmap, tt, xyzRange); 
[e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);       

JB = J*B;
[Np, Nt] = size(JB);
Np = length(weightIm(mask));
W = spdiag(weightIm(mask),0,Np,Np);
sqW = spdiag(sqrt(weightIm(mask)),0,Np,Np);
Wv = weightIm(mask);
Wvmat = repmat(Wv,1,Nt);
sqWvmat = repmat(sqrt(Wv),1,Nt);
dkB = 2*real(1i*JB'*W*e);

grad_norm0 = norm(dkB);
grad_norm = grad_norm0;
b = dpls(As,d,beta,mask);
b0 = b; 
cost0 = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
nrmse0 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));

betavec = sqrt(beta)*ones(size(b));
%% projected GD

beta_linSearch = 0.5; % parameter for backtracking line search
dobreak = 0;
niterK = 5;   % number of k and b update iterations before resetting momentum
exetime = 0;
coeff = coeff0(:);
quadOptions = optimset('Algorithm', 'interior-point-convex','Display','off');
[coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -double(coeff), Call, 0.99*bound,[],[],[],[],[],quadOptions);
k = basis*reshape(coeff,[length(coeff)/3,3]);
[A,mask,xx,yy,zz,nx,ny,nz] = formA(k, sens, mask, fieldmap, tt, xyzRange); 
Ain = A;
kin = k;

%%%%% the following few lines record the error right after B-spline fitting
[e, As, ~] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
b = dpls(As,d,beta,mask);
cost1 = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
nrmse1 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
k1 = k; 
b1 = b; 

printf('iter | cost | nrmse | nrmse_ov')
t = 1.0;  % step size
coeff_old = coeff;

quadOptions = optimset('Algorithm', 'interior-point-convex','Display','off','MaxIter',20);

fistaMode = 0;   % 0: FISTA, 1: MFISTA

deltak=1;
k_allIter{2} = k;
iter = 1;
nrmseold = 1;
costold = 1000;
deltacost = 1000;
deltanrmse = 1;
exetime_start = toc;
while iter < (niter+1) & deltacost > 0 %deltanrmse > 5e-3 % deltacost > 0.2 % deltak > 0.05
%for iter = 1:niter
   % optimize k
   
	% wipe out momentum from previous iteration
	% check out MFISTA?
	momentum = 1; %min(1 + 4/(10^3)*(iter-1)^3, 2)
	coeff_y = coeff;   % momentum term ('y' in the FISTA paper)

	for iterK = 1:niterK
		
		%startit=toc;

		% Step 4.1 in FISTA paper: Minimize quadratic Taylor approximation at current value of coeff_y ("z" in barista paper)
		
		% get the gradient and J'J matrix for k
		kv = B*coeff_y;
		k = reshape(kv,size(k));
		A = formA(k, sens, mask, fieldmap, tt, xyzRange); 
		%A = updateA(Ain,k-kin,mask,xx,yy,zz,nx,ny,nz);
		[e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
		%fprintf(1,'kb_Jacobian_input(): %f sec\n', toc-start);
		JB = J*B;
		WJB = Wvmat.*JB;  % this is about 25% faster than W*JB.
		sqWJB = sqWvmat.*JB;
		%dkB = 2*(2*pi)*real(1i*JB'*W*e);   % gradient
		%HB = 2*(2*pi)^2*real(JB'*W*JB);    % Hessian
		dkB = 2*(2*pi)*real(1i*WJB'*e);   % gradient
		HB = 2*(2*pi)^2*real(sqWJB'*sqWJB);    % Hessian. This is about 2x faster than (JB'*W*JB), see http://www.walkingrandomly.com/?p=4912

		% Do the quadratic constrained optimization (step 4.1 in FISTA paper)
		[deltacoeff,FVAL,EXITFLAG,OUTPUT] = quadprog(HB, dkB, Call, 0.99*(bound-Call*coeff_y),[],[],[],[],[],quadOptions);
		%fprintf(1,'quadprog(): %f sec\n', toc-start);
		coeff_old = coeff;               % x_{k-1} in 4.3 in FISTA paper. Solution estimate from previous iteration.
		switch fistaMode
			case 0
				coeff = coeff + t*deltacoeff;      % x_k in 4.1 in FISTA paper, or x_{k=1} in barista paper
			case 1
				coeff_z = coeff + t*deltacoeff;    % z_k above 5.2 in MFISTA paper.
		end
      
		% evaluate the cost function
		kv = B*coeff;
		k = reshape(kv,size(k));
		%A = updateA(Ain,k-kin,mask,xx,yy,zz,nx,ny,nz);
		A = formA(k, sens, mask, fieldmap, tt, xyzRange); 
		[e, As] = eval_func(A, b, d, sens, mask);
		cost = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
		
		% rough b update
   	%b = qpwls_pcg1_hao(b, As, W, d(mask), diag_sp(betavec),'niter',2);
    b = bupdate(b, As, W, d, mask, betavec,2);

		% update momentum (steps 4.2 and 4.3 in FISTA paper) if momentum direction is not too far off current direction. See barista paper.
		momentum_old = momentum; % nesterov momentum
		momentum = 1; %(1 + sqrt(1 + 4*momentum_old^2))/2;
		coeff_y_old =  coeff_y;
		coeff_y = coeff + (momentum_old - 1)/momentum*(coeff - coeff_old);
		dcy = (coeff(:)-coeff_y_old(:));
		dc  = (coeff(:)-coeff_old(:));
		kappa = norm(dcy)*norm(dc);
		if ( dcy'*dc / kappa) < cos(10/180*pi) 
			% wipe out momentum
			coeff_y = coeff;
			momentum = 1;
			%fprintf(1,'momentum reset to 1 at iterK = %d \n', iterK);
		end
	
		switch fistaMode
			case 0
			case 1
				% evaluate cost function at previous iterate, and at coeff_z (z_k in MFISTA paper)
				coeff
				% evaluate the cost function
				kv = B*coeff;
				k = reshape(kv,size(k));
				%A = updateA(Ain,k-kin,mask,xx,yy,zz,nx,ny,nz);
				A = formA(k, sens, mask, fieldmap, tt, xyzRange); 
				[e, As] = eval_func(A, b, d, sens, mask);
				cost = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
		end

		%fprintf(1,'one iteration: %f sec\n', toc-startit);
   end % end of niterK operations
   
   % b update 
   % b2 = dpls(A,d,beta,mask);
   %b = qpwls_pcg1_hao(b, As, W, d(mask), diag_sp(betavec),'niter',10);
   b = bupdate(b, As, W, d, mask, betavec, 10);

   % relax phase within iv_roi (actually doesn't seem to help suppress OV excitation so no need for this)
   %d_st = As*b;   % small-tip approximation of excitation pattern
   %d(mask) = d(mask)./exp(1i*angle(d(mask))).*exp(1i*angle(d_st));  % it doesn't matter that we change the target phase outside the IV as well

   % record the residual error after update b.
   e = d(mask) - As*b; 

   % record the residual error in OV
   projIm = embed(As*b,mask);   % [32 32 12]
   %save projIm projIm
   e_ov = projIm(ov_roi);   % OV excitation 
   nrmse_ov(iter) = norm(e_ov)/sqrt(numel(e_ov))/max(abs(d(mask)));
   %norm(e_ov)     %/max(abs(d(mask)))
   
   %%
   deltak = norm(k_allIter{iter+1}(:)-k(:),1)/size(k,1);  % mean change in k per sample [cycles/fov I think]
   k_allIter{iter+2} = k;
   b_allIter{iter+2} = b; 
   costRec(iter) = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
   deltacost = abs(costRec(iter)-costold)/costold*100; %  percent change
   costold = costRec(iter);
   nrmse(iter) = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
   deltanrmse = abs(nrmse(iter)-nrmseold);
   nrmseold = nrmse(iter);
   exetime(iter) = toc-exetime_start;
   printf('%d | %.5g | %.5g | %.5g | %.3g \n', iter, costRec(iter), nrmse(iter), nrmse_ov(iter));
   
   if grad_norm/grad_norm0 < 0.001
      break;
   end
   if dobreak == 1;
      break;
   end
  	iter = iter+1; 
end
k_allIter{1} = k0;
b_allIter{1} = b0;
k_allIter{2} = k1;
b_allIter{2} = b1;
costRec = [cost0 cost1 costRec]; 
nrmse = [nrmse0 nrmse1 nrmse]; 
exetime = [0 0 exetime]; 
save mfista_exetime costRec nrmse nrmse_ov exetime


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
