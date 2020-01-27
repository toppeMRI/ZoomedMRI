function [b1,gx,gy,gz,projIm,M_pr,K_rv] = ktContRF_greedy2(d,fieldmap,sens, mask,tipangle,n_pe,ishard,xyzRange,tol,redo)
% find ktpoints using greedy approach (daehyun's method) and then connect them using fastest gradient approach. Works well.

if ~isvar('tol') % parameter for cubic spline smoothing
    tol = 0.2;
end
if ~isvar('redo') % parameter for cubic spline smoothing
    redo = 1;
end
doGD = 0; % do local minimization using Levenberg-Marquardt algorithm.
redo = 0; % refind all phase encoding locations
optAfter = 0;
% load single_coil_tx_sensitivity;
% sens = ones(size(fieldmap));  % in the spins trajectory experiment, I used all ones as sens
% load sensitivities;
sens = permute(sens, [2,3,4,1]); % ncoil is the first dimension expect choose PE code
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
% zfov = 20;

dx = 64; % number of pixles
dy = 64;
dz = length(xyzRange.z);
% dz = 5;
% pulse parameters for spoke pulse design
% pulse_params.kz_area = 6/slice_thickness;            %-4 ~ 4 cycle/cm : 4/slice_thickness or more recommended for a good slice profile
pulse_params.tip_angle = tipangle;         %10 degree
% pulse_params.slice_thickness = slice_thickness;  %
pulse_params.dt = 4*(10^-6);         %4usec
pulse_params.slew_rate = 15000;      %g/cm/sec, max slew rate
pulse_params.g_max = 4;              % 4g, maximum gradient amp.
pulse_params.max_pe_jump = 15/xfov;

%% find kt points using greedy approach
% if 1
if redo == 1
    [nrmse,kx,ky,kz] = determine_pe_field_greedy3_3d(d,sens,mask,n_pe,type,ncand,0*fieldmap,xfov,yfov,zfov,pulse_params,0,ishard);
    % [nrmse,kx,ky,kz] = determine_pe_field_inverseA_3d(d,sens,mask,n_pe,type,ncand,fieldmap,xfov,yfov,zfov,pulse_params,0,ishard);
    [kx, ky, kz]
end
% end

%% gradient descent
load kxyztmp;
load weighting;
k = [kx,ky,kz];
k_old = k;

if doGD
    tic
    %    load kxyztmp;
    %    [k, cost] = PSD_ktraj(k, d, sens, mask, weightIm,beta);
    [k, costQN] = Newton_LM_ktraj_full(k_old, d, sens, mask, weightIm,beta);
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
    k_old = [k_old; zeros(1,3)];
    n_pe = n_pe + 1;
    doGD = 1;
    save kxyztmp;
    toc
end

%% Connect those points using tsp algorithm
% if 1 % for testing
if redo == 1
    addpath C:\Users\Hao\Dropbox\MRI_research\pulseDesign\SPINS\ktpoints\tsp_ga;
    addpath ~/Dropbox/MRI_research/pulseDesign/SPINS/ktpoints/tsp_ga;
    addpath .\tsp_ga;
    % tsp_ga;  demo
    xyz = k;
    popSize = 60;
    numIter = 5e3;
    showProg = 1;
    showResult = 1;
    
    a = meshgrid(1:n_pe);
    dmat = reshape(sqrt(sum((xyz(a,:)-xyz(a',:)).^2,2)),n_pe,n_pe);
    [optRoute,minDist] = tsp_ga(xyz,dmat,popSize,numIter,showProg,showResult);
    kstart = find(optRoute == n_pe);
    optRoute = [optRoute(kstart+1:end), optRoute(1:kstart)];
    kx_order = k(optRoute,1);
    ky_order = k(optRoute,2);
    kz_order = k(optRoute,3);
    kx_old_order = k_old(optRoute,1);
    ky_old_order = k_old(optRoute,2);
    kz_old_order = k_old(optRoute,3);
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
    save kxyztmp; % save to avoid run the previous code everytime in debugging.
end
% end



%% Find the shortest time gradient given the path
addpath(genpath('C:\Users\Hao\Dropbox\MRI_research\pulseDesign\SPINS\ktpoints\minTimeGradient\'));
addpath(genpath('~/Dropbox/MRI_research/pulseDesign/SPINS/ktpoints/minTimeGradient/'));
% addpath(genpath('.\minTimeGradient'));
tolTmp = tol;
if doGD == 0
    load kxyztmp;
    k = [kx_old_order,ky_old_order,kz_old_order];
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

ks = [kx/xfov,ky/yfov,kz/zfov];
[C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(ks,0, 0, 0, 4, 15, 4e-3);     % Rotationally variant solution

%% gradient descent
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
    [ks, cost] = Newton_LM_ktraj_full(k(1:end-endk,:), d, sens, mask, weightIm,beta);
    [ks, cost] = ADMM_ktraj(k(1:end-endk,:), d, sens, mask, weightIm,beta,xfov,yfov,zfov);
    ksall(:,1) = [ks(:,1); k(end-endk+1:end,1)]/xfov;
    ksall(:,2) = [ks(:,2); k(end-endk+1:end,2)]/yfov;
    ksall(:,3) = [ks(:,3); k(end-endk+1:end,3)]/zfov;
    [C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(ksall,0, 0, 0, 4, 15, 4e-3);     % Rotationally variant solution
    save kfinal k ks
end


%% gradient descent
if optAfter % hard to choose a proper step size
    ndt = 10;
    k = K_rv(1:ndt:end,:);
    [ks, cost] = fmincon_ktraj(k, d, sens, mask, weightIm,beta);
    [C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(ks,0, 0, 0, 4, 15, 4e-3);     % Rotationally variant solution
    save kfinal k ks
end

%% plot for the grant application
figure,plot3(kx_old_order,ky_old_order,kz_old_order,'+');
hold on,plot3(kx_order,ky_order,kz_order,'r^');
hold on,plot3(kx_order,ky_order,kz_order,'g--');
hold on,plot3(C_rv(:,1)*xfov, C_rv(:,2)*yfov, C_rv(:,3)*zfov,'k-');
grid on ; xlabel('kx'),ylabel('ky'),zlabel('kz');
%% Generate RF
d = d.*sin(tipangle/180*pi);
% [b1,gx,gy,gz] = ktContRF(d,fieldmap,mask,tipangle,n_pe,ishard,xyzRange)

sens = permute(sens,[4,1,2,3]);
load sensitivities;
% G_rv = zeros(size(G_rv));
% K_rv = zeros(size(K_rv));
[b1,gx,gy,gz,projIm] = cont_gen_rf(d,fieldmap,mask,sens,xyzRange,K_rv,G_rv);
%figure,im(projIm); title('small tip approximation image');


%% plot simulations and pulse diagram

% sens = permute(sens,[4,1,2,3]);
T1 = 1000; T2 = 100; dt = 4e-6;
M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt,fieldmap,mask,T1/1000,T2/1000);
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

