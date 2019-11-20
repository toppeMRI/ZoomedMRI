function [ko_pcm, optInfos] = BsplInteriorPoint(phm, d, ki_pcm, varargin)
% Interior Point Method w/ 2nd-order B-spline on kTraj
% The output k-traj of this method is not necessarity balanced, but ends at 0
% Ref. Sun et al. DOI: 10.1109/TMI.2015.2478880
% This function needs MIRT
%INPUTS
% - phm
%   .sMap
%   .w
%   .fov
%   .b0Map
% - d
% - ki_pcm (nk, ndim), cycle/cm
%OPTIONAL
% - alpha    (1,), dflt 0.01, parameter for backtracking line search
% - beta     (1,), dflt 0.5, parameter for backtracking line search
% - eta      (1,), dflt 8, coeff for $\|rf\|_2^2$ cost
% - ncycle
% - niter_o
% - niter_i  (1,), # of 
% - niter_rf (1,), # of cg-steps computing rf
% - niter_nt (1,), # of Newton steps fixing log-barrier cost
% - dt       (1,), Sec, dwell time
% - gam      (1,), Hz/G, gyro freq
% - gMax     (1,), G/cm, max gradient
% - sMax     (1,), G/cm/Sec, max slew rate of gradient
% - isFatrix

%% init
% all parameters empirically chosen by Sun
arg = envMR('get_s','gMax','sMax','dt','gam');
arg.isFatrix = false;
[arg.alpha, arg.beta, arg.eta, arg.ncycle] = deal(0.01, 0.5, 0.006, 100);
[arg.niter_o, arg.niter_i, arg.niter_rf, arg.niter_nt] = deal(20, 2, 80, 2);

arg = attrParser(arg, varargin);

[alpha, beta, eta, ncycle] = deal(arg.alpha, arg.beta, arg.eta, arg.ncycle);
[niter_o,  niter_i]  = deal(arg.niter_o,  arg.niter_i);
[niter_rf, niter_nt] = deal(arg.niter_rf, arg.niter_nt);
[gMax, sMax, dt, gam] = deal(arg.gMax, arg.sMax, arg.dt, arg.gam);

FT_c = {'dt',dt, 'gam',gam, 'isFatrix',arg.isFatrix};

[sMap, w, fov, b0Map, ofst] = getattrs(phm, {'sMap','w','fov','b0Map','ofst'});
m = ~~w;
[dm, wm] = deal(d(m), w(m));
k_pcm = ki_pcm;

%% parameterize by 2nd-order bspline, and form the peak loc bounds
[nk, ndim] = size(k_pcm);
[B, dB1, dB2] = kTx.basis.Bspline2ndOrder_regular(nk, ncycle, ndim, dt, gam);
[dB_k, bnd]=deal([dB1;dB2],[gMax*ones(size(dB1,1),1);sMax*ones(size(dB2,1),1)]);
B_k = kron(speye(3), spCompatible(B)); % _k: kron'ed

coeff = B\k_pcm; % cycle/cm/B
coeff_v = spCompatible(coeff(:));

%% get feasible initial for IPM
% project 2nd-order bspline parameterized traj to the feasible set
coeff_v = quadprog(eye(numel(coeff_v)), -coeff_v, dB_k, bnd);
coeff = reshape(coeff_v, size(B,2), []);
k_pcm = B*coeff;

%% Prep
kpcm2k_fn  = @(kpcm)bsxfun(@times, kpcm, fov);
A_fn       = @(k)form_pTxst(form_FTst(k,m,b0Map, 'ofst',ofst,FT_c{:}), sMap);
%{
rf_fn      = @(A, rf0)qpwls_pcg1(rf0,A,Gdiag(wm),dm, eta*size(A,2)...
                                 ,'niter',niter_rf);
%}
rf_fn      = @(A, dummy)(A'*A + eta*size(A,2)*speye(size(A,2)))\A'*dm;
em_fn      = @(A, rf)(dm - A*rf);
% cost_fn    = @(em,rf)real(em'*(wm.*em) + eta*size(rf,1)*(rf'*rf));
cost_fn    = @(em,rf)real(em'*(wm.*em)); % + eta*size(rf,1)*(rf'*rf));
nrmse_fn   = @(em)norm(em)/sqrt(numel(em))/max(abs(dm));
barrier_fn = @(coeff_v)(bnd - dB_k*coeff_v);

% after parameterize by 2nd-order bspline
k  = kpcm2k_fn(k_pcm);
A  = A_fn(k);
rf = rf_fn(A, []);
% rf = (A'*A + speye(size(A,2)))\A'*dm; % init first `rf` in a robust way
em = em_fn(A, rf);
cost1  = cost_fn(em, rf);
nrmse1 = nrmse_fn(em);

[k1, rf1] = deal(k, rf);

%% Interior point
mu_int = 8; % multiply to update the weighting t_int

dobreak = 0;

tic;

for iter = 1:niter_o
  t_int = 80;
  % optimize k
  for iInter = 1:niter_i
    % printf('iteration of interior point, %d\n',iInter)
    t_int = t_int*mu_int; % trade off between the data fitting and log-barier
    for iNewton = 1:niter_nt
      % get the gradient and Hessian matrix for k
      [~, ~, ~, J] = formJ_ls(A, rf, m);
      JB = J*B_k;
      em = em_fn(A, rf);

      dkB = 2*(2*pi)   * real(1i*JB'*(wm.*em));
      HB  = 2*(2*pi)^2 * real(JB'*(bsxfun(@times, wm, JB)));
      % gradient and hessian of log barier function
      % Below adds in the log barrier part
      barrier  = barrier_fn(coeff_v);
      d_phi    = 1./barrier;
      grad_phi = dB_k'*d_phi;
      hess_phi = dB_k'*(bsxfun(@times, d_phi.^2, dB_k));
      
      ipdkB = dkB*t_int + grad_phi;
      ipHB  = HB *t_int + hess_phi;
      dir   = spCompatible(ipHB\ipdkB); % descending direction
      
      phi     = -sum(log(barrier)); % log(prod(barrier)) fails due to precision
      ipCost  = cost_fn(em, rf)*t_int + phi;
      ipCost0 = ipCost;
      
      % backtrack line search: reduce step-size if the ipCost doesn't decrease
      % enough or the new solution is infeasible.
      t = 1;
      coeff_v0 = coeff_v;
      dcost = dir'*ipdkB;
      while (ipCost >= ipCost0-alpha*t*dcost) || any( barrier_fn(coeff_v)<0 )
        % update coeffs
        coeff_v = coeff_v0 - t*dir;
        
        k_pcm = reshape(B_k*coeff_v, size(k)); % update k_pcm
        k  = kpcm2k_fn(k_pcm);
        A  = A_fn(k); % update k traj (and A) only, fix rf
        em = em_fn(A, rf);
        
        phi = -sum(log(barrier_fn(coeff_v)));
        ipCost = cost_fn(em,rf)*t_int + phi;
        t = beta * t; % reduce step size
        
        % break if step size is too small
        dobreak = all(abs(t * dir) < 1e-5); % break if step size is too small
        if dobreak, warning('step size is too small'); break; end
      end
      disp([string(num2str(cost_fn(em,rf)))])
      
      % printf('ip_iter, %.1f %.5f, %.5f\n',t_int, t, cost_fn(em,rf));
      if(dobreak == 1), break; end
    end % end of newton iterations for one iterior point on the central path
  end % end of all interior points
  
  rf = rf_fn(A, rf); % update b with newest k, cg faster than direct inverse
  em = em_fn(A, rf); % record the residule error after update rf.
  
  % k_allIter{iter+1}  = k;
  % rf_allIter{iter+1} = rf;
  [cost(iter), nrmse(iter)] = deal(cost_fn(em, rf), nrmse_fn(em));
  if ~mod(iter,5),printf('%d | %.5g | %.5g \n',iter,cost(iter),nrmse(iter)); end
  
  if dobreak == 1, break; end
end

% k_allIter{1}  = k1;
% rf_allIter{1} = rf1;
cost = [cost1 cost];
nrmse = [nrmse1 nrmse];

ko_pcm = k_pcm;

optInfos.costRec = cost;
optInfos.nrmse   = nrmse;
% optInfos.k_allIter = k_allIter;
% optInfos.rf_allIter = rf_allIter;
end

