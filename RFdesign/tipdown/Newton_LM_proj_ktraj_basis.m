%% Newton method full update (precondition gradient descent with backtracking)
% full update to k, alternate update b
% Levenberg-Marquardt method; add lambda*I to hessian and adjust lambda

function [k,costA, nrmse] = Newton_LM_proj_ktraj_basis(k, d, sens,mask,weightIm,beta,niter,fieldmap,tt, xyzRange, ndt, xfov, yfov, zfov);

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


%% projected LM
printf('iter | cost | nrmse | lam | norm(g)')
dobreak = 0;
lam = 1e-10;
coeff = coeff0(:);
alp0 = 1;
alp = alp0;

quadOptions = optimset('Algorithm', 'active-set'); 
[coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -coeff, Call, 0.99*bound,[],[],[],[],[],quadOptions); 
k = basis*reshape(coeff,[length(coeff)/3,3]); 

for iter = 1:niter
    % get the Hessian matrix for k
    kv = k(:);
    [H] = getHessian(k,b,d,A,e, sens,mask,weightIm,xyzRange);
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
        
        % projection
        [coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -coeff, Call, 0.99*bound,[],[],[],[],[],quadOptions); 

        
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


function [H] = getHessian(k,b,d,A,e,sens,mask,weightIm,xyzRange)
%% construct matrices

x = xyzRange.x/(xyzRange.x(2)-xyzRange.x(1));
y = xyzRange.y/(xyzRange.y(2)-xyzRange.y(1));
z = xyzRange.z/(xyzRange.z(2)-xyzRange.z(1));
[xx,yy,zz]=ndgrid(x,y,z);

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
