% Interior point algorithm
function [k, costRec, nrmse, exetime, k_allIter, b_allIter] = InteriorPoint_ktraj_basis(k, d, sens,mask,weightIm,beta,niter,fieldmap,tt, xyzRange, ndt, xfov, yfov, zfov)
%
% k = [Nt 3] = [kx, ky, kz]

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
t_real = 4e-3*[0:length(k)-1]; 
basis = bspline_1d_synth(coefs,tis,'ending','zero','order',2); figure(2),plot(t_real, basis);
% basis = single(basis);
% basis = eye(length(k));
B = kron(eye(3,3),basis);
B = sparse(B); 
coeff0 = basis\k;
figure(2),plot(k),hold on,plot(basis*coeff0,'--');
k0 = k; 

[Nt, Nl] = size(basis);

%% construct difference matrix and bound
nk = size(k,1);
C1 = Cdiff1(nk,'order',1);  % [Nt Nt]
C2 = Cdiff1(nk,'order',2);
C1_basis = C1 * basis;   % [Nt Nl]. This is "DH" in the paper.
C2_basis = C2 * basis;   % [Nt Nl]

[Y, maxIdx] = max(C1_basis);   % [1 Nl]
[Y, minIdx] = min(C1_basis);
peakIdxC1 = union(maxIdx,minIdx);   % These are the rows that "P_1" in the paper picks out.
peakIdxC2 = peakIdxC1 + round((peakIdxC1(2) - peakIdxC1(1))/3);  % Corresponds to "P_2" in the paper.
nps = length(peakIdxC1); % number of peaks
C1_basis_peak = C1_basis(peakIdxC1,:);   % This is "PDH" in the paper (Eq. (6)).
C2_basis_peak = C2_basis(peakIdxC2,:);
% Call = [C2_basis_peak]; % relax the gradient constraint
Call = [C1_basis_peak; C2_basis_peak];
Call = double([Call; -Call]);  % Each block in Eq. (6) in the paper, except ordered/organized a bit differently.
%Call = [Call,zeros(size(Call)),zeros(size(Call)); zeros(size(Call)), Call, zeros(size(Call)); zeros(size(Call)), zeros(size(Call)), Call];
Call = kron(eye(3),Call);  % "U" in the paper, Eq. (6)

gamma = 4257.6; % Hz/Gauss
Gmax = 3.9; % Gauss/cm
Smax = 14.9e3; % Gauss/cm/s
dt = 4e-6; %s
bound1 = Gmax*gamma*dt*ndt; %*[ones(nps,1)*xfov,ones(nps,1)*yfov,ones(nps,1)*zfov]; % cycle/fov
bound2 = Smax*gamma*(dt)^2*ndt^2; %*[ones(nps,1)*xfov,ones(nps,1)*yfov,ones(nps,1)*zfov];
boundX = [bound1*ones(nps,1)*xfov; bound2*ones(nps,1)*xfov; bound1*ones(nps,1)*xfov; bound2*ones(nps,1)*xfov];
boundY = [bound1*ones(nps,1)*yfov; bound2*ones(nps,1)*yfov; bound1*ones(nps,1)*yfov; bound2*ones(nps,1)*yfov];
boundZ = [bound1*ones(nps,1)*zfov; bound2*ones(nps,1)*zfov; bound1*ones(nps,1)*zfov; bound2*ones(nps,1)*zfov];

% % if relax the gradient constraint
% boundX = [bound2*ones(nps,1)*xfov; bound2*ones(nps,1)*xfov;];
% boundY = [bound2*ones(nps,1)*yfov; bound2*ones(nps,1)*yfov;];
% boundZ = [bound2*ones(nps,1)*zfov; bound2*ones(nps,1)*zfov;];
 
bound = [boundX; boundY; boundZ]; % positive and negative;   This is Eq. (8) in the paper.


%% get initial
b = 0;
% [dkx, dky, dkz, db, e, As, J] = kb_grad_withB0(k, b, d, sens, mask,weightIm,beta,fieldmap,tt,xyzRange);
% [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
A = formA(k0, sens, mask, fieldmap, tt, xyzRange);
[e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);

% J is size [Np 3Nt], B is size [3Nt 3L], where L = # basis functions
JB = J*B;   % size [Np 3L]. This is 'J' = [Jx,Jy,Jz] in the paper, with Jx as in Eq. (10) etc.

Np = length(weightIm(mask));
W = spdiag(weightIm(mask),0,Np,Np);
dkB = 2*real(1i*JB'*W*e);

grad_norm0 = norm(dkB);
grad_norm = grad_norm0;
b = dpls(As,d,beta,mask);
b0 = b; 
cost0 = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
nrmse0 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));

betavec = sqrt(beta)*ones(size(b));

%% Interior point
% gradient of log barrier:

nInter = 2;  % number of steps on the central path
nNewton = 2; % number of newton iterations for each step
mu_int = 2; % multiply to update the weighting t_int 

printf('iter | cost | nrmse | lam | norm(g)')
dobreak = 0;
lam = 1e-10;
coeff = coeff0(:);   % 'coeff' is 'c' in the paper, i.e., the basis coefficients that we wish to solve for.

alp = 1;

quadOptions = optimset('Algorithm', 'interior-point-convex','Display','off');
[coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -double(coeff), Call, double(0.99*bound),[],[],[],[],[],quadOptions);
k = basis*reshape(coeff,[length(coeff)/3,3]);
A = formA(k, sens, mask, fieldmap, tt, xyzRange); 

%%%%% the following few lines record the error right after B-spline fitting
[e, As, ~] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
b = dpls(As,d,beta,mask);
cost1 = norm(sqrt(W)*e)^2 + beta*norm(b)^2;
nrmse1 = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
k1 = k; 
b1 = b; 

exetime = 0;
tic;
alpha = 0.01; % parameter for backtracking line search
beta_linSearch = 0.5; % parameter for backtracking line search


for iter = 1:niter

	% b is updated in this loop 

   t_int = 20;

   % Optimize coeff (k), holding current estimate of b fixed.
   for iInter = 1:nInter
		
		% t_int is updated in this loop

      printf('iteration of interior point, %d\n',iInter)
      t_int = t_int*mu_int; % parameter controls the trade off between the data fitting and log-barier

      for iNewton = 1:nNewton

			% coeff is updated in this loop (holding t_int and b fixed)

         % Get the gradient and Hessian matrix, using the most recent update of coeff (and hence A).
         
         [e, As, J] = kb_Jacobian_inputA(A, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
         % [e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
         JB = J*B;   % Eq. (10) in the paper.
         dkB = 2*(2*pi)*real(1i*JB'*W*e);  % Eq. (11) in the paper.
         HB = 2*(2*pi)^2*real(JB'*W*JB);   % Eq. (14) in the paper.
         
         d_phi = 1./(bound - Call*coeff); % 'w' in Eqs. (16-17) in the paper. Gradient and hessian of log barrier function. 
         grad_phi = Call'*d_phi;  % Eq. (16)
         hess_phi = Call'*spdiag(d_phi.^2, 0, Nt, Nt)*Call;  % Eq. (17)
         
         dkBAll = dkB*t_int + grad_phi;
         HBAll = HB*t_int + hess_phi;

			% descending direction. This is 'delta' in Algorithm 2, step 4 in the paper.
			dire = (HBAll)\dkBAll; 
         
         coeff_old = coeff;
         
         phi = -sum(log(bound - Call*coeff));
         cost = (norm(sqrt(W)*e)^2 + beta*norm(b)^2)*t_int + phi;
         cost_old = cost;
         
         % backtrack line search: reduce step if the cost doesn't decrease
         % enough or the new solution (coeff) is infeasible. 
         t = 1;
         while (cost == cost_old - alpha*t*(dire'*dkBAll)) || (cost > cost_old - alpha*t*(dire'*dkBAll)...
               || (min(bound - Call*coeff) < 0) )  % the last condition checks for feasibility

            % update coeff (and k-space trajectory)
            coeff = coeff_old -  t * dire;
            kv = B*coeff;  % 'B' here is 'H' in the paper (basis function matrix).
            k = reshape(kv,size(k));
            
            % evaluate the cost function, Eq. (15)
            A = formA(k, sens, mask, fieldmap, tt, xyzRange);
            [e, As] = eval_func(A, b, d, sens, mask);
            phi = -sum(log(bound - Call*coeff));
            cost = (norm(sqrt(W)*e)^2 + beta*norm(b)^2)*t_int + phi;  % Eq. (15) in the paper.

				% reduce step size (to be used next, if current coeff is infeasible or doesn't reduce cost)
            t = beta_linSearch * t; 
            
            % break if step size is too small
            maxgo = max(abs(t * dire));
            if maxgo < 0.00001
               dobreak = 1;
               warning('step size is too small');
               break;
            end
         end % end of while

			% We've exited while loop, i.e., we've found a coeff we're happy with (in this Newton iteration).
			% Note that A has been updated according to the current coeff.
         
         printf('NW,%.1f %.5f, %.5f, %d, %d, %d\n',t_int, t, (norm(sqrt(W)*e)^2 + beta*norm(b)^2), cost == cost_old - alpha*t*(dire'*dkBAll),(cost > cost_old - alpha*t*(dire'*dkBAll)), (min(bound - Call*coeff) < 0));
         if(dobreak == 1)
            break;
         end
      end % end of newton iterations for one iterior point on the central path
   end % end of all interior points
   
	% Update b
   %    b2 = dpls(A,d,beta,mask); 
   b = qpwls_pcg1_hao(b, As, W, d(mask), diag_sp(betavec),'niter',20);
   e = d(mask) - As*b; % record the residual error after update b.
   
   %%
   k_allIter{iter+2} = k; % recorded only for plot since the simulation resolution is 64x64
   b_allIter{iter+2} = b; 
   costRec(iter) = norm(sqrt(W)*e)^2 + beta*norm(b)^2; 
   nrmse(iter) = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
   exetime(iter) = toc;
   printf('%d | %.5g | %.5g | %.5g | %.3g \n', iter, costRec(iter), nrmse(iter), lam, grad_norm)
   lam = lam/10;
   
   %    if grad_norm/grad_norm0 < 0.05
   %       break;
   %    end
   if dobreak == 1;
      break;
   end
end
k_allIter{1} = k0;
b_allIter{1} = b0;
k_allIter{2} = k1;
b_allIter{2} = b1;
costRec = [cost0 cost1 costRec]; 
nrmse = [nrmse0 nrmse1 nrmse]; 
exetime = [0 0 exetime]; 
save ip_exetime exetime costRec nrmse
end


