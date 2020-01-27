function [b1,gx,gy,gz,projIm,M_pr,K_rv] = ktContRF_greedy2(d,fieldmap,fieldmap_simu,sens, mask,tipangle,n_pe,ishard,xyzRange,tol,extype,isBssfp,reply)
% find ktpoints using greedy approach (daehyun's method) and then connect them using fastest gradient approach. Works well.

useSingle = 0;
if useSingle == 1
   d = single(d);
   fieldmap = single(fieldmap);
   sens = single(sens);
   tipangle = single(tipangle);
   xyzRange.x = single(xyzRange.x);
   xyzRange.y = single(xyzRange.y);
   xyzRange.z = single(xyzRange.z);
end

if ~isvar('tol') % parameter for cubic spline smoothing
   tol = 0.2;
end
if ~isvar('redo') % refind all phase encoding locations
   redo = 1;
end
if ~isvar('isBssfp') % make excitation balanced
   isBssfp = 0;
end
doGD = 0; % do local minimization using Levenberg-Marquardt algorithm.
redo = 1; % refind all phase encoding locations
noOrder = 0; % doesn find the shortest path visiting order (original KT-points);
% load single_coil_tx_sensitivity;
% sens = ones(size(fieldmap));  % in the spins trajectory experiment, I used all ones as sens
% load sensitivities;
sens = permute(sens, [2,3,4,1]); % ncoil is the first dimension except choosing PE code; it is not changed back later because I am working on single coil now.
mask = logical(mask);
beta = 8; % regularizar on b

ncand = 10; %# of candidate in KT points puse design
type = 2;

xfov = 24;
yfov = 24;
if length(xyzRange.z) > 1
   zfov = xyzRange.z(end) - xyzRange.z(1) + xyzRange.z(2) - xyzRange.z(1);
else
   zfov = 0;
end

% pulse parameters for spoke pulse design
% pulse_params.kz_area = 6/slice_thickness;            %-4 ~ 4 cycle/cm : 4/slice_thickness or more recommended for a good slice profile
pulse_params.tip_angle = tipangle;         %10 degree
% pulse_params.slice_thickness = slice_thickness;  %
pulse_params.dt = 4*(10^-6);         %4usec
pulse_params.slew_rate = 15000;      %g/cm/sec, max slew rate
pulse_params.g_max = 4;              % 4g, maximum gradient amp.
pulse_params.max_pe_jump = 15/xfov;

maxSlewInMinTimeGradient = 15; % g/cm/ms  The fastest gradient method using optimal control can has small violation of constraints.

%% find kt points using greedy approach
if 0

if redo == 1
   [nrmse,kx,ky,kz] = determine_pe_field_greedy3_3d(d,sens,mask,n_pe,type,ncand,0*fieldmap,xfov,yfov,zfov,pulse_params,0,ishard);
   %    [nrmse,kx,ky,kz] = determine_pe_field_inverseA_3d(d,mask, n_pe);
   [kx, ky, kz]
end
% end


%% Connect those points using tsp algorithm
% if 1 % for testing
if redo == 1
   if noOrder == 0
      addpath C:\Users\Hao\Dropbox\MRI_research\pulseDesign\SPINS\ktpoints\tsp_ga;
      addpath ~/Dropbox/MRI_research/pulseDesign/SPINS/ktpoints/tsp_ga;
      addpath .\tsp_ga;
      % tsp_ga;  demo
      xyz = [kx,ky,kz];
      popSize = 20;
      numIter = 5e3;
      showProg = 1;
      showResult = 1;
      a = meshgrid(1:n_pe);
      dmat = reshape(sqrt(sum((xyz(a,:)-xyz(a',:)).^2,2)),n_pe,n_pe);
      [optRoute,minDist] = tsp_ga(xyz,dmat,popSize,numIter,showProg,showResult);
      kstart = find(optRoute == n_pe);
      optRoute = [optRoute(kstart+1:end), optRoute(1:kstart)];
      kx_order = xyz(optRoute,1);
      ky_order = xyz(optRoute,2);
      kz_order = xyz(optRoute,3);
      %     kx_old_order = k_old(optRoute,1);
      %     ky_old_order = k_old(optRoute,2);
      %     kz_old_order = k_old(optRoute,3);
      knorm = sum([kx_order,ky_order,kz_order].^2,2);
      nonZeroK_ind = find(knorm~=0);
      
      if knorm(nonZeroK_ind(end)) > knorm(nonZeroK_ind(1))
         kx_order((nonZeroK_ind(1):nonZeroK_ind(end))) = flipud(kx_order(nonZeroK_ind(1):nonZeroK_ind(end)));
         ky_order((nonZeroK_ind(1):nonZeroK_ind(end))) = flipud(ky_order(nonZeroK_ind(1):nonZeroK_ind(end)));
         kz_order((nonZeroK_ind(1):nonZeroK_ind(end))) = flipud(kz_order(nonZeroK_ind(1):nonZeroK_ind(end)));
      end
      kx_order(1:max(nonZeroK_ind(1)-1,1)) = [];
      ky_order(1:max(nonZeroK_ind(1)-1,1)) = [];
      kz_order(1:max(nonZeroK_ind(1)-1,1)) = [];
      k = [kx_order,ky_order,kz_order];
   else
      kx_order = kx;
      ky_order = ky;
      kz_order = kz;
   end
   save kxyztmp; % save to avoid run the previous code everytime in debugging.
end
% end
%% gradient descent
% unit for k here: cycles/fov
if doGD
   tic
   load kxyztmp;
   load weighting;
   k = [kx_order,ky_order,kz_order];
   %    [k, cost] = PSD_ktraj(k, d, sens, mask, weightIm,beta);
   k_old = k;
   [k, costGD] = Newton_LM_ktraj_full(k_old, d, sens, mask, weightIm,beta, 20, xyzRange);
   %    [k, costQN] = QN_ktraj(k_old, d, sens, mask, weightIm,beta);
   %    [k, cost] = Newton_ktraj(k_old, d, sens, mask, weightIm,beta);
   %    [k, costFull] = Newton_ktraj_full(k_old, d, sens, mask, weightIm,beta);
   %    [k, costkb] = Newton_ktraj_kb(k_old, d, sens, mask, weightIm,beta);
   %    [k, cost] = SD_ktraj(k, d, sens, mask, weightIm,beta);
   %    figure,plot(cost); hold on, plot(costFull,'r--');
   %    hold on, plot(costQN,'g-o')
   %    hold on, plot(costkb,'k.-');
   %    legend('alternate (seperate P)','alternate (full P)', 'full update of kb (seperate P)', 'QN');
   %    kx_order = k(:,1);
   %    ky_order = k(:,2);
   %    kz_order = k(:,3);
   k = [k; zeros(1,3)];
   doGD = 1;
   save kxyztmp;
   toc
   figure,plot(log10(costGD));
end


%% Find the shortest time gradient given the path
addpath(genpath('C:\Users\Hao\Dropbox\MRI_research\pulseDesign\SPINS\ktpoints\minTimeGradient\'));
addpath(genpath('K:\Dropbox\Dropbox\MRI_research\pulseDesign\SPINS\ktpoints\minTimeGradient\'));
addpath(genpath('~/Dropbox/MRI_research/pulseDesign/SPINS/ktpoints/minTimeGradient/'));
% addpath(genpath('.\minTimeGradient'));
tolTmp = tol;
if doGD == 0
   load kxyztmp;
   k = [kx_order,ky_order,kz_order];
else
   load kxyztmp;
end
tol = tolTmp;
% kx = kx_order;
% ky = ky_order;
% kz = kz_order;
kx = k(:,1);
ky = k(:,2);
kz = k(:,3);

w = ones(size(kx)); w([1 end]) = 100; %smooth the kt points
[sp1, kx_sp] = spaps([1:length(kx)],kx,tol,w,3);
[sp2, ky_sp] = spaps([1:length(kx)],ky,tol,w,3);
[sp3, kz_sp] = spaps([1:length(kx)],kz,tol,w,3);
figure,plot(kx),hold on,plot(kx_sp,'r');
kx = kx_sp';
ky = ky_sp';
kz = kz_sp';

ks = [kx/xfov,ky/yfov,kz/zfov];   % only ks, K_rv are in cycle/cm because of the convention in fastest gradient algorithm
[C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(ks,0, 0, 0, 4, maxSlewInMinTimeGradient, 4e-3);     % Rotationally variant solution

% % The fastest gradient method using optimal control can has small violation
% of constraints. Run this several times helps.
% [C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(K_rv,0, 0, 0, 4, maxSlewInMinTimeGradient, 4e-3);
% [C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(K_rv,0, 0, 0, 4, maxSlewInMinTimeGradient, 4e-3);

%% project to feasible K
maxSlewInMinTimeGradient = 15; % since I save and load variables in this code, so the value may not be 15 before this line
if max(abs(S_rv(:))) > maxSlewInMinTimeGradient
   K_rv_old = K_rv;
   K_rv = projK(K_rv,maxSlewInMinTimeGradient);
   G_rv = k2g_hao(K_rv);
   dt = 4e-6; 
   S_rv = diff(G_rv,1,1)/(dt*1000); S_rv = [S_rv; [0,0,0]];
end

%% optimize directly on the final k-traj (very slow)
%
if 0 % hard to choose a proper step size
   K_rv_old = K_rv;
   load weighting;
   % k in this part is in 1/cm
   kx = K_rv(:,1)*xfov; % for convergence of GD algorithm
   ky = K_rv(:,2)*yfov;
   kz = K_rv(:,3)*zfov;
   k = [kx,ky,kz];
   k = k(1:4:end,:);
   endk = 0; % keep the end part of the traj fixed.
   %       [ks, cost] = Newton_LM_ktraj_full(k(1:end-endk,:), d, sens, mask, weightIm,beta);
   %    [ks, cost] = Newton_ktraj(k(1:end-endk,:), d, sens, mask, weightIm,beta);
   [ks, cost] = ADMM_ktraj(k(1:end-endk,:), d, sens, mask, weightIm,beta,xfov,yfov,zfov);
   ksall(:,1) = [ks(:,1); k(end-endk+1:end,1)]/xfov;
   ksall(:,2) = [ks(:,2); k(end-endk+1:end,2)]/yfov;
   ksall(:,3) = [ks(:,3); k(end-endk+1:end,3)]/zfov;
   [C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(ksall,0, 0, 0, 4, maxSlewInMinTimeGradient, 4e-3);     % Rotationally variant solution
   save kfinal k ks
end



%%
%%
%% optimize directly on the final k-traj (using basis function)
%%
%%

%     K_rv = zeros(size(K_rv)); % test for prephasing problem
%     G_rv = zeros(size(G_rv)); % test for prephasing problem

% load k_sos; % stack of spiral initial pulse
% load k_spins; % stack of spiral initial pulse...............
% load k_ktpoints_prephasing; 
% load k_ktpoints_rFOV; 
% load k_extendedKT_prephasing2; 
% load k_extendedKT_rFOV; 
% load test; 

% K_rv = k3d;
% G_rv = g3d;

% K_rv = K_rv(end-981:end,:); % 3.9 ms % length of k_extendedKT_rFOV
% G_rv = G_rv(end-981:end,:); 

save test K_rv G_rv; 
else
	load test
end

dt = 4e-6;
S_rv = diff(G_rv,1,1)/(dt*1000); S_rv = [S_rv; [0,0,0]];


fprintf(1,'\n');
if ~exist('reply','var')
	reply = input('Run interior point optimization? (y/n) ', 's');
end
fprintf(1,'\n');
if strcmp(reply,'y')   % To generate IV pulse for first rFOV paper, don't run interior-point optimization (JFN)
%if 1
   ndt = 1;
   load weighting;
   K_rv_old = K_rv;
   k = K_rv(1:ndt:end,:);
   dt = 4e-6;
   tt = [0:length(k)-1]*ndt*dt-length(k)*ndt*dt;
   endk = 0;
   %     [ks, cost] = Newton_LM_ktraj_full(k, d, sens, mask, weightIm,beta,niter);
   kx = k(:,1)*xfov; % for convergence of GD algorithm cycle/fov
   ky = k(:,2)*yfov;
   kz = k(:,3)*zfov;
   k = [kx,ky,kz];
	save k k
   niter = 20;
   ks = k;
   
   method = 3;
   %    [ks, cost] = Newton_LM_ktraj_full(k, d, sens, mask, weightIm,beta,niter);
   %    [ks, cost] = Newton_LM_ktraj_slewRate(k, d, sens, mask, weightIm,beta,niter,xfov,yfov,zfov,4e-3*ndt);
   %    [ks, cost] = Newton_LM_ktraj_penalized(k, d, sens, mask, weightIm,beta,niter);
   
   %    [ks, cost, nrmse] = Newton_LM_ktraj_basis(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange); % work
   %    [ks, cost, nrmse] = Newton_LM_proj_ktraj_basis(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
   %     [ks, cost, nrmse, exetime] = Newton_LM_proj_ktraj_basis_lookUpHessian(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
   %      [ks, cost, nrmse, exetime] = Newton_LM_proj_ktraj_basis_Kanzow(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
   
   %     [ks, cost, nrmse, exetime] = InteriorPoint_ktraj_basis(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
   
   %    [ks, cost, nrmse] = BBGM_proj_ktraj_basis(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % not work with projection
   %    [ks, cost, nrmse] = GD_proj_ktraj_basis(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
   %    [ks, cost, nrmse] = fmincon_ktraj_basis(k, d, sens, mask, weightIm, beta,niter,fieldmap, ndt, tt, xfov, yfov, zfov,xyzRange); % work
   %    [ks, cost] = fmincon_ktraj(k, d, sens, mask, weightIm, beta,niter,fieldmap, ndt, tt, xfov, yfov, zfov); % has gradient/slew rate vialation, check later.
   %    figure,plot(cost);title('cost function vs iteration');
   switch method
      case 1
         [ks, cost, nrmse, exetime, k_allIter, b_allIter] = GD_proj_ktraj_basis(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
      case 2
         %          [ks, cost, nrmse, exetime] = Newton_LM_proj_ktraj_basis_lookUpHessian(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
         [ks, cost, nrmse, exetime, k_allIter, b_allIter] = Newton_LM_proj_ktraj_basis_Kanzow(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov); % work
      case 3
			% Perform Iterated Local Search
			clear cost
			save fov xfov yfov zfov   % needed in projK.m
			mxs = 15;  % max slew rate. G/cm/msec

			ii=1;
			if 0
	      	[kin, costRec, nrmse, exetime, k_allIter, b_allIter] = ...
					InteriorPoint_ktraj_basis(k,d,sens,mask,weightIm,beta,niter,fieldmap,tt,xyzRange,ndt,xfov,yfov,zfov); 
				cost(ii,:) = costRec;  ii = ii+1;
				k = kin;
				%for nCycle = [25 26 38 33 49] % for nCycle=50 in Interpoint_Ktraj_basis.m, this sequence produces best solution.
				for nCycle = [] % for nCycle=50 in Interpoint_Ktraj_basis.m, this sequence produces best solution.
					kin = projK(kin,mxs,nCycle);  
	      		[kin, costRec, nrmse, exetime, k_allIter, b_allIter] = ...
						InteriorPoint_ktraj_basis(kin,d,sens,mask,weightIm,beta,niter,fieldmap,tt,xyzRange,ndt,xfov,yfov,zfov); 
					cost(ii,:) = costRec;  ii = ii+1;
				end
				save kin
			else
				load kin
				k = kin;
			end

			% exhaustive search using 'kin' as starting point
			nCycle = [31 33 29 27];  % for nCycle=100 in Interpoint_ktraj_basis.m
			nCycle = [30 25 20 15 10]; % for nCycle=50 in Interpoint_Ktraj_basis.m 25 is best.
			nCycle = [35:-1:15]; % for nCycle=50 in Interpoint_Ktraj_basis.m 25 is best.
			nCycle = [49:-1:20]; % 
			nCycle = []; % 
			for it = 1:length(nCycle)
				nCycle(it)

				% Perturb k-space to (hopefully) a different basin of attraction.
				k = projK(kin,mxs,nCycle(it));    

				% Perform local search
	         [k, costRec, nrmse, exetime, k_allIter, b_allIter] = ...
					InteriorPoint_ktraj_basis(k, d, sens, mask, weightIm,beta,niter,fieldmap,tt,xyzRange,ndt, xfov, yfov, zfov);
				cost(ii,:) = costRec;  ii = ii+1;
			end
			ks = k;
			save cost cost
      case 4
         [ks, cost, nrmse, exetime, k_allIter, b_allIter] = fmincon_ktraj_basis(k, d, sens, mask, weightIm, beta,niter,fieldmap, ndt, tt, xfov, yfov, zfov,xyzRange); % work
   end
   exetime(end)
	cost = cost(end);
   kstmp(:,1) = [ks(:,1); k(end-endk+1:end,1)]/xfov;
   kstmp(:,2) = [ks(:,2); k(end-endk+1:end,2)]/yfov;
   kstmp(:,3) = [ks(:,3); k(end-endk+1:end,3)]/zfov;
   %    figure,plot([0:length(nrmse)-1], nrmse); title('nrmse vs iteration');

   % for plot
   save(['convergeSpeed' num2str(method)], 'exetime', 'nrmse', 'cost', 'k_allIter', 'b_allIter', 'd', 'sens', 'mask','fieldmap', 'tt', 'xyzRange');
   figure,plot(exetime, nrmse); title('nrmse vs time'); xlabel('time (sec)'); ylabel('NRMSE');
   figure,plot(exetime, cost); title('cost vs time'); xlabel('time (sec)'); ylabel('a.u.');
   
   ks = [kstmp];
   
   if extype == 2 || extype == 4 || extype == 3 % include 3 here for debug, delete later
      ksLength = size(ks,1);
      K_rv = interp1(1:ksLength, ks, 1:1/ndt:ksLength);
      G_rv = k2g_hao(K_rv);
      S_rv = diff(G_rv,1,1)/(dt*1000); S_rv = [S_rv; [0,0,0]];
   else
      [C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(ks,0, 0, 0, 4, maxSlewInMinTimeGradient, 4e-3);     % Rotationally variant solution
   end
   save kfinal K_rv ks K_rv_old G_rv
   
   t = [0:size(K_rv,1)-1]*dt*1000; %ms for plot
   figure,
   subplot(3,1,1),plot(t, K_rv(:,:)); hold on; plot(t, K_rv_old(:,:),'--');legend('kx','ky','kz'); ylabel('cycle/cm'); title('k-space trajectory');
   subplot(3,1,2),plot(t, G_rv); hold on; plot(t, 4*ones(size(G_rv)),'k--'); hold on; plot(t, -4*ones(size(G_rv)),'k--'); ylim([-4.2, 4.2]); legend('gx','gy','gz','limit'); ylabel('gauss/cm'); title('gradient waveform');
   subplot(3,1,3),plot(t, S_rv); hold on; plot(t, 15*ones(size(S_rv)),'k--'); hold on; plot(t, -15*ones(size(S_rv)),'k--'); ylim([-20, 20]); legend('sx','sy','sz','limit'); ylabel('gauss/cm/mm'); title('slew rate');
   xlabel('time (ms)');
   
end


%% project to feasible K
maxSlewInMinTimeGradient = 15; % since I save and load variables in this code, so the value may not be 15 before this line
if max(abs(S_rv(:))) > maxSlewInMinTimeGradient
   K_rv_old = K_rv;
   K_rv = projK(K_rv,maxSlewInMinTimeGradient);
   G_rv = k2g_hao(K_rv);
   dt = 4e-6; 
   S_rv = diff(G_rv,1,1)/(dt*1000); S_rv = [S_rv; [0,0,0]];
   
   t = [0:size(K_rv,1)-1]*dt*1000; %ms for plot
   figure,
   subplot(3,1,1),plot(t, K_rv_old(:,:)); hold on; plot(t, K_rv(:,:),'--');legend('kx','ky','kz','k-proj'); ylabel('cycle/cm'); title('before/after projection: k-space trajectory');
   subplot(3,1,2),plot(t, G_rv); hold on; plot(t, 4*ones(size(G_rv)),'r--');legend('gx','gy','gz','limit'); ylabel('gauss/cm'); title('gradient waveform');
   subplot(3,1,3),plot(t, S_rv); hold on; plot(t, 15*ones(size(S_rv)),'r--'); legend('sx','sy','sz','limit'); ylabel('gauss/cm/mm'); title('slew rate');
   xlabel('time (sec)');
end

%% add a ramp to the beginning of gradient so slew rate is not violated when adding gradient crusher
G_rv_begin = G_rv(1,:); 
ramp_nt = ceil(max(abs(G_rv_begin))./maxSlewInMinTimeGradient/4e-3)+1; 
ramp_x = linspace(0, G_rv_begin(1), ramp_nt); 
ramp_y = linspace(0, G_rv_begin(2), ramp_nt); 
ramp_z = linspace(0, G_rv_begin(3), ramp_nt);
G_rv = [[ramp_x', ramp_y', ramp_z']; G_rv];
G_rv = [zeros(1,3); G_rv]; 
K_rv = g2k_hao(G_rv);
S_rv = diff(G_rv,1,1)/(dt*1000);




%% make it balanced
if isBssfp == 1
   sumG = sum(G_rv,1);
   ncycles = 0;
   opslthick = 0.5; % this value doesn't matter...
   gxRephaser = makecrusher_hao(ncycles, opslthick, sumG(1)*4e-6);
   gyRephaser = makecrusher_hao(ncycles, opslthick, sumG(2)*4e-6);
   gzRephaser = makecrusher_hao(ncycles, opslthick, sumG(3)*4e-6);
   maxGLength = max([length(gxRephaser), length(gyRephaser), length(gzRephaser)]);
   gxRephaser = [gxRephaser; zeros(maxGLength - length(gxRephaser),1)];
   gyRephaser = [gyRephaser; zeros(maxGLength - length(gyRephaser),1)];
   gzRephaser = [gzRephaser; zeros(maxGLength - length(gzRephaser),1)];
   G_rv = [[gxRephaser, gyRephaser, gzRephaser]; G_rv];
   K_rv = g2k_hao(G_rv);
   S_rv = diff(G_rv,1,1)/(dt*1000);
end

%% plot for the grant application
if 0
   figure,plot3(ks(:,1),ks(:,2),ks(:,3),'o');
   hold on,plot3(K_rv(:,1), K_rv(:,2), K_rv(:,3),'k-');
   grid on ; xlabel('kx'),ylabel('ky'),zlabel('kz');
end

%% Generate RF
% d = d.*sin(tipangle/180*pi);
sens = permute(sens,[4,1,2,3]); % put coil in the first dimension
load sensitivities;
load weighting;
% [b1,gx,gy,gz,projIm] = cont_gen_rf(d,fieldmap,mask,sens,xyzRange,K_rv,G_rv);
%figure,im(projIm); title('small tip approximation image');

% use direct method here:
dt = 4e-6; tt = [0:length(K_rv)-1]*dt-length(K_rv)*dt;
A = formA([K_rv(:,1)*xfov, K_rv(:,2)*yfov, K_rv(:,3)*zfov], permute(sens,[2,3,4,1]), mask, fieldmap, tt, xyzRange);

% combine sens and A into As
ncoils = size(sens,1);
if ncoils ==1
   As = A;
else
   Atmp = [];
   for ic = 1:ncoils
      sensi = sens(ic,:,:,:);
      Atmp = [Atmp; spdiag(sensi(mask))*A];
   end
   As = reshape(Atmp, [sum(mask(:)),ncoils*length(kx)]);
end

betavec = sqrt(beta)*ones(size(As,2),1);
W = diag_sp(weightIm(mask));
b1 = qpwls_pcg1_hao(zeros(size(betavec)), As, W, d(mask), diag_sp(betavec),'niter',100);
% b1 = dpls(As,d,beta,mask);
% addpath minMax
% b1 = min_split_ADMM(As, d(mask)); 
g = k2g_hao(K_rv);
gx = g(:,1);
gy = g(:,2);
gz = g(:,3);
projIm = embed(As*b1,mask);

% scale with flip angle:
b1 = b1.*sin(tipangle/180*pi);
projIm = projIm.*sin(tipangle/180*pi);


%% plot simulations and pulse diagram

% sens = permute(sens,[4,1,2,3]);
T1 = 1000; T2 = 100; dt = 4e-6;
Nt = length(gx);
M_pr = parallel_blochCim(0,double(reshape(b1,[Nt,ncoils])),double(gx),double(gy),double(gz),double(sens),double(xyzRange.x),double(xyzRange.y),double(xyzRange.z),dt,double(fieldmap_simu),mask,T1/1000,T2/1000);

% figure,im(abs(M_pr));colormap default;
tt = 1:length(b1);
dt = 4e-3;
tt = tt*dt;
% figure,im('row',2,cat(3,projIm, M_pr));colormap default;

% figure,im(angle(M_pr./d)/pi*180);colormap default;
% figure,subplot(4,1,1),plot(tt, abs(b1(:,1))); ylabel('gauss');
% subplot(4,1,2),plot(tt,abs(gx)); ylabel('g/cm');
% subplot(4,1,3),plot(tt,abs(gy));
% subplot(4,1,4),plot(tt,abs(gz));

