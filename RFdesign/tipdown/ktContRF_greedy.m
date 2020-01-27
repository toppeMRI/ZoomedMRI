function [b1,gx,gy,gz] = ktContRF_greedy(d,fieldmap,mask,tipangle,n_pe,ishard,xyzRange)
% find ktpoints using greedy approach (daehyun's method) and then connect them using fastest gradient approach. Works well. 

tic;

% load single_coil_tx_sensitivity;
sens = ones(size(fieldmap));  % in the spins trajectory experiment, I used all ones as sens
mask = logical(mask); 


ncand = 10;
type = 2;

ncoils = 1; !!!
lambda = 5; % penalty parameter for spoke pulse
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


[nrmse,kx,ky,kz] = determine_pe_field_greedy3_3d(d,sens,mask,n_pe,type,ncand,fieldmap,xfov,yfov,zfov,pulse_params,0,ishard);
% [nrmse,kx,ky,kz] = determine_pe_field_inverseA_3d(d,sens,mask,n_pe,type,ncand,fieldmap,xfov,yfov,zfov,pulse_params,0,ishard);
ks = [kx/xfov,ky/yfov,kz/zfov]; 

addpath(genpath('minTimeGradient/')); 
[C_rv, time_rv, G_rv, S_rv, K_rv] = minTimeGradient(ks,1, 0, 0, 4, 15, 4e-3);     % Rotationally variant solution


d = d.*sin(tipangle/180*pi); 
% [b1,gx,gy,gz] = ktContRF(d,fieldmap,mask,tipangle,n_pe,ishard,xyzRange)
sens = ones([1, size(d)]); 
[b1,gx,gy,gz] = cont_gen_rf(d,fieldmap,mask,sens,xyzRange,C_rv,G_rv); 

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
