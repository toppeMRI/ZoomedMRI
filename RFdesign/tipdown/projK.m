%% project the K traj to feasible region
function k = projK(k,Smax,nCycle)

kin = k;

if ~exist('nCycle', 'var')
	nCycle = 100;
end

gamma = 4257.6; % Hz/Gauss
Gmax = 3.9; % Gauss/cm
Smax = Smax*1000; % Gauss/cm/s
dt = 4e-6; %s

if 0
   ndt = 1;
   nk = size(k,1);
   C1 = Cdiff1(nk,'order',1);
   C2 = Cdiff1(nk,'order',2);
   Call = double(full([C1; C2]));
   Call = [Call; -Call];
   Call = [Call,zeros(size(Call)),zeros(size(Call)); zeros(size(Call)), Call, zeros(size(Call)); zeros(size(Call)), zeros(size(Call)), Call];
   
   bound1 = Gmax*gamma*dt*ndt; %*[ones(nk,1)*xfov,ones(nk,1)*yfov,ones(nk,1)*zfov]; % cycle/fov
   bound2 = Smax*gamma*(dt)^2*ndt^2; %*[ones(nk,1)*xfov,ones(nk,1)*yfov,ones(nk,1)*zfov];
   boundX = [bound1*ones(nk,1); bound2*ones(nk,1); bound1*ones(nk,1); bound2*ones(nk,1)];
   boundY = [bound1*ones(nk,1); bound2*ones(nk,1); bound1*ones(nk,1); bound2*ones(nk,1)];
   boundZ = [bound1*ones(nk,1); bound2*ones(nk,1); bound1*ones(nk,1); bound2*ones(nk,1)];
   bound = [boundX; boundY; boundZ];
   
   % project coeff to feasible region first
   quadOptions = optimset('Algorithm', 'active-set');
   [kv,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(k(:))), -k(:), Call, 0.99*bound,[],[],[],[],[],quadOptions);
   k = reshape(kv, k);
else
   ndt = 1;
   %nCycle = 200;
   coefs = eye(nCycle);
   iEnd = nCycle + 1;
   ti = linspace(0,iEnd,size(k,1))';
   tis = repmat(ti, [1,nCycle]);
   basis = bspline_1d_synth(coefs,tis,'ending','zero','order',2); %figure(2),plot(ti, basis);
   basis = double(basis);
   B = kron(eye(3,3),basis);
   
   coeff = basis\k;
   %figure(2),plot(k),hold on,plot(basis*coeff,'--'); pause(1);
   kold = k;
   k = basis*coeff;
   coeff = coeff(:);
   
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
   
	load fov
   bound1 = Gmax*gamma*dt*ndt; %*[ones(nps,1)*xfov,ones(nps,1)*yfov,ones(nps,1)*zfov]; % cycle/fov
   bound2 = Smax*gamma*(dt)^2*ndt^2; %*[ones(nps,1)*xfov,ones(nps,1)*yfov,ones(nps,1)*zfov];
   boundX = [bound1*ones(nps,1); bound2*ones(nps,1); bound1*ones(nps,1); bound2*ones(nps,1)]*xfov;
   boundY = [bound1*ones(nps,1); bound2*ones(nps,1); bound1*ones(nps,1); bound2*ones(nps,1)]*yfov;
   boundZ = [bound1*ones(nps,1); bound2*ones(nps,1); bound1*ones(nps,1); bound2*ones(nps,1)]*zfov;
   bound = [boundX; boundY; boundZ]; % positive and negative
   
   %quadOptions = optimset('Algorithm', 'active-set');
   quadOptions = optimset('Algorithm', 'interior-point-convex');
   [coeff,FVAL,EXITFLAG,OUTPUT] = quadprog(eye(length(coeff)), -coeff, Call, 0.99*bound,[],[],[],[],[],quadOptions);
   k = basis*reshape(coeff,[length(coeff)/3,3]); 
end

