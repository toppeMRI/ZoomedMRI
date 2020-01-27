function [b1,gx,gy,gz] = ktContRF(d,fieldmap,mask,tipangle,n_pe,ishard,xyzRange)
% find ktpoints and they construct continuous wave from those points 
% those ktpoints can be find from greedy approach or inverse system matrix
% approach or inverse DFT approach. Try inverse system matrix first. 
save tmp
load tmp
slew_max = 15000;      %g/cm/sec, max slew rate
g_max = 4;              % 4g, maximum gradient amp.

% load single_coil_tx_sensitivity;
sens = ones(size(fieldmap));  % in the spins trajectory experiment, I used all ones as sens
mask = logical(mask); 

%% find kt points using inverse system matrix approach. 
[nrmse,kx,ky,kz] = determine_pe_field_inverseA_3d(d,mask,n_pe); 


%% Connect those points using tsp algorithm
addpath ./tsp_ga; 
% tsp_ga;  demo

xyz = [kx,ky,kz];
popSize = 60;
numIter = 2e4; 
showProg = 1;
showResult = 1;
a = meshgrid(1:n_pe);
dmat = reshape(sqrt(sum((xyz(a,:)-xyz(a',:)).^2,2)),n_pe,n_pe);
[optRoute,minDist] = tsp_ga(xyz,dmat,popSize,numIter,showProg,showResult);
save tmp; 
%% do some processing to the TSP trajectory
load tmp; 
kstart = find(optRoute == n_pe); 
optRoute = [optRoute(kstart+1:end), optRoute(1:kstart)]; 
kx_orded = kx(optRoute); 
ky_orded = ky(optRoute); 
kz_orded = kz(optRoute); 
kx_orded(1:8) = []; 
ky_orded(1:8) = []; 
kz_orded(1:8) = []; 
% figure,subplot(3,1,1), plot(kx_orded); 
% subplot(3,1,2), plot(ky_orded); 
% subplot(3,1,3), plot(kz_orded); 
figure,plot3(kx_orded,ky_orded,kz_orded); title('k-traj after smooth'); 

%% smooth k-trajectory
Nsmo = round(n_pe/5);
Ntrail = 3; % # of points padding to the end of k-traj to reach to zero
x_orded = 1:length(kx_orded); 
x_smo = linspace(1,x_orded(end),Nsmo);
%kx_smo = fit(x_smo',kx_orded, 'smoothingspline'); 
poly_smo = polyfit(x_orded', kx_orded, round(n_pe/15));
kx_smo = polyval(poly_smo,x_smo);
poly_smo = polyfit(x_orded', ky_orded, round(n_pe/15));
ky_smo = polyval(poly_smo,x_smo); 
poly_smo = polyfit(x_orded', kz_orded, round(n_pe/15));
kz_smo = polyval(poly_smo,x_smo);
kx_smo = [kx_smo,linspace(kx_smo(end),0,Ntrail)]';
ky_smo = [ky_smo,linspace(ky_smo(end),0,Ntrail)]';
kz_smo = [kz_smo,linspace(kz_smo(end),0,Ntrail)]';
x_smo = [x_smo, [x_smo(end) + 1: x_smo(end) + Ntrail]];  
k_smo = [kx_smo, ky_smo, kz_smo]; 
figure,
subplot(3,1,1), plot(kx_orded); hold on, plot(x_smo, kx_smo,'r'); 
subplot(3,1,2), plot(ky_orded); hold on, plot(x_smo, ky_smo,'r'); 
subplot(3,1,3), plot(kz_orded); hold on, plot(x_smo, kz_smo,'r'); 
figure,plot3(kx_smo,ky_smo,kz_smo); title('k-traj after smooth'); 

% %% generate the fastest gradient based on Sana Vaziri, and Michael Lustig
% addpath(genpath('minTimeGradient\')); 
% [C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(k_smo,1, 0, 0, 4, 15, 4e-3);     % Rotationally variant solution

%% interpolation k-trajectory to desired pulse length
Trf = 1e-3; %s
dt = 4e-6; %s 
Nrf = Trf/dt; 
so = 1:length(kx_smo); 
si = linspace(1,so(end),Nrf); 
kx_smo = interp1(so,kx_smo,si,'spline'); 
ky_smo = interp1(so,ky_smo,si,'spline'); 
kz_smo = interp1(so,kz_smo,si,'spline'); 
% figure,
% subplot(3,1,1), plot(kxi); title('kx'); 
% subplot(3,1,2), plot(kyi); title('ky'); 
% subplot(3,1,3), plot(kzi); title('kz'); title('k-traj after interpolation'); 
% kxi = smooth(kxi)'; 
% kyi = smooth(kyi)'; 
% kzi = smooth(kzi)'; 
k_smo = [kx_smo', ky_smo', kz_smo']; 
% figure,plot3(kxi,kyi,kzi); title('k-traj after interpolation'); 

fovx = xyzRange.x(end) - xyzRange.x(1) + xyzRange.x(2) - xyzRange.x(1); 
fovy = xyzRange.y(end) - xyzRange.y(1) + xyzRange.y(2) - xyzRange.y(1);  
if length(xyzRange.z) > 1
    fovz = xyzRange.z(end) - xyzRange.z(1) + xyzRange.z(2) - xyzRange.z(1);
else
    fovz = 0.5; % 2d imaging, this number doesn't matter as long as it is nonzero
end

save tmp; 
load tmp; 

%% generate gradient using the normal way
Ng = length(kx_smo); 
gam=26751;                   %rad/sec/g
gambar = gam/2/pi;           %Hz/g
dt = 4e-6; 
kG = -gambar*dt*triu(ones([Ng,Ng])); 
kGi = inv(kG); 
ak = [kx_smo'/fovx, ky_smo'/fovy, kz_smo'/fovz]; %Hz/cm 
G = []; 
G = kGi * ak;   %G/cm
Gx = G(:,1); 
Gy = G(:,2); 
Gz = G(:,3); 
S = []; 
S(:,1) = gradient(G(:,1),dt); %Hz/cm/sec
S(:,2) = gradient(G(:,2),dt); %Hz/cm/sec
S(:,3) = gradient(G(:,3),dt); %Hz/cm/sec
figure,subplot(3,1,1), plot(G(:,1)); ylabel('Gx'); 
hold on; plot(4*ones(size(G,1)),'r'); 
hold on; plot(-4*ones(size(G,1)),'r'); 
subplot(3,1,2), plot(G(:,2)); ylabel('Gy'); 
hold on; plot(4*ones(size(G,1)),'r'); 
hold on; plot(-4*ones(size(G,1)),'r');
subplot(3,1,3), plot(G(:,3)); ylabel('Gz'); 
hold on; plot(4*ones(size(G,1)),'r'); 
hold on; plot(-4*ones(size(G,1)),'r');
legend('gradient (G/cm)', 'gradient limit (4 G/cm)'); 

figure,subplot(3,1,1), plot(S(:,1)); ylabel('Sx'); 
hold on; plot(15000*ones(size(S,1)),'r'); 
hold on; plot(-15000*ones(size(S,1)),'r'); 
subplot(3,1,2), plot(S(:,2)); ylabel('Sy'); 
hold on; plot(15000*ones(size(S,1)),'r'); 
hold on; plot(-15000*ones(size(S,1)),'r');
subplot(3,1,3), plot(S(:,3)); ylabel('Sz'); 
hold on; plot(15000*ones(size(S,1)),'r'); 
hold on; plot(-15000*ones(size(S,1)),'r');
legend('slew rate (G/cm)', 'slew rate limit (15000 G/sec/cm)'); 

save tmp; 
load tmp; 


%% smooth lead-out of gradient
nfit = 8; % # of smaples to fit inside gradient
npad = 5; % pad leading samples with a few zeros (for PCG)
ntrail = 3; % pad zeros to the end of gradient wave
t_lead = max(G(1,:))/slew_max; % minimum # samples for lead out
nlead = ceil(t_lead/dt); 
Glead = []; 
for ii = 1:3
%   Glead(:,ii) =  [zeros([npad 1]);linspace(0,G(1,ii),nlead)'];
    Glead(:,ii) =  zeros(npad+nlead,1);
end
G = cat(1,Glead,G); 
S = []; 
S(:,1) = gradient(G(:,1),dt); %Hz/cm/sec
S(:,2) = gradient(G(:,2),dt); %Hz/cm/sec
S(:,3) = gradient(G(:,3),dt); %Hz/cm/sec

figure,subplot(3,1,1), plot(S(:,1)); ylabel('Sx'); 
hold on; plot(15000*ones(size(S,1)),'r'); 
hold on; plot(-15000*ones(size(S,1)),'r'); 
subplot(3,1,2), plot(S(:,2)); ylabel('Sy'); 
hold on; plot(15000*ones(size(S,1)),'r'); 
hold on; plot(-15000*ones(size(S,1)),'r');
subplot(3,1,3), plot(S(:,3)); ylabel('Sz'); 
hold on; plot(15000*ones(size(S,1)),'r'); 
hold on; plot(-15000*ones(size(S,1)),'r');
legend('slew rate (G/cm)', 'slew rate limit (15000 G/sec/cm)'); 

%% iterative smooth
Ng = length(G(:,1)); 
while abs(max(S(:))) > 15000
    G(:,1) = smooth(G(:,1)); 
    G(:,2) = smooth(G(:,2)); 
    G(:,3) = smooth(G(:,3));
    G(:,1) = G(:,1) - sum(G(:,1))/Ng; 
    G(:,2) = G(:,2) - sum(G(:,2))/Ng; 
    G(:,3) = G(:,3) - sum(G(:,3))/Ng; 
    S = zeros(size(G)); 
    S(:,1) = gradient(G(:,1),dt); %Hz/cm/sec
    S(:,2) = gradient(G(:,2),dt); %Hz/cm/sec
    S(:,3) = gradient(G(:,3),dt); %Hz/cm/sec
end
ntrail = 4; 
for ii = 1:3
%   Glead(:,ii) =  [zeros([npad 1]);linspace(0,G(1,ii),nlead)'];
Gtrail(:,ii) = linspace(G(end,ii),0,ntrail); 
end
G = [G; Gtrail]; 
Nrf = size(G,1); 
kG = -gambar*dt*triu(ones([Nrf,Nrf])); 
k_smo = kG*G; 
kx = k_smo(:,1); 
ky = k_smo(:,2); 
kz = k_smo(:,3); 
figure,plot3(kx,ky,kz); 
save tmp;
%% Generate RF waves
load tmp; 
d = d.*sin(tipangle/180*pi); 
% [b1,gx,gy,gz] = ktContRF(d,fieldmap,mask,tipangle,n_pe,ishard,xyzRange)
sens = ones([1, size(d)]); 
[b1,gx,gy,gz] = cont_gen_rf(d,fieldmap,mask,sens,xyzRange,k_smo,G); 

%% plot simulations and pulse diagram 

% sens = permute(sens,[4,1,2,3]);  
T1 = 1000; T2 = 100; dt = 4e-6;
M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt,fieldmap,mask,T1/1000,T2/1000);
% figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default;
figure,im(abs(M_pr));colormap default;
tt = 1:length(b1);
dt = 4e-3; 
tt = tt*dt; 
figure,subplot(4,1,1),plot(tt, abs(b1(:,1))); ylabel('gauss'); 
subplot(4,1,2),plot(tt,abs(gx)); ylabel('g/cm'); 
subplot(4,1,3),plot(tt,abs(gy)); 
subplot(4,1,4),plot(tt,abs(gz)); ylim([-5,5]); 


end




