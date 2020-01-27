function [b1,gx,gy,gz] = ktpointsRF(d,fieldmap,mask,tipangle,n_pe,ishard,xyzRange)
% ktpoints using greedy approach (daehyun's method). Works well. 
% !! set fieldmap to zero when select phase encoding locations
% load single_coil_tx_sensitivity;
sens = ones(size(fieldmap));  % in the spins trajectory experiment, I used all ones as sens
mask = logical(mask); 


ncand = 10;
type = 2;

ncoils = 1; !!!
lambda = 5; % penalty parameter for spoke pulse
xfov = 24;
yfov = 24;
zfov = xyzRange.z(end) - xyzRange.z(1) + xyzRange.z(2) - xyzRange.z(1);
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

[nrmse,kx,ky,kz] = determine_pe_field_greedy3_3d(d,sens,mask,n_pe,type,ncand,fieldmap,xfov,yfov,zfov,pulse_params,0,ishard,0);
% [nrmse,kx,ky,kz] = determine_pe_field_inverseA_3d(d,sens,mask,n_pe,type,ncand,fieldmap,xfov,yfov,zfov,pulse_params,0,ishard);
ks = [kx,ky,kz];
[x_est,projection,pp_nrmse] = determine_weights_field2_3d(xfov,yfov,zfov,kx,ky,kz,d,mask,sens,fieldmap,pulse_params,lambda,0,ishard);

[b1,gx,gy,gz] = create_pulse_3d(x_est,n_pe,ncoils,kx,ky,kz,pulse_params,dx,dy,dz,xfov,yfov,zfov,ishard);
b1 = b1*1e4; % convert to gauss
disp(sprintf('the net pulse length %f msec',size(b1,1)*pulse_params.dt*1000));

% x_range = ([0:1:dx-1]-floor(dx/2))*xfov/dx;
% y_range = ([0:1:dy-1]-floor(dy/2))*yfov/dy;
% z_range = [0];
% 

g3d = [gx, gy, gz]; 
k3d = g2k_hao(g3d); 
save k_ktpoints_rFOV g3d k3d
%% plot simulations and pulse diagram 
sens = permute(sens,[4,1,2,3]);  
T1 = 1000; T2 = 100; dt = 4e-6;
M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt,fieldmap,mask,T1/1000,T2/1000);
figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default;
figure,im(abs(M_pr));colormap default;

tt = 1:length(b1);
dt = 4e-3; 
tt = tt*dt; 
figure,subplot(4,1,1),plot(tt, abs(b1(:,1))); ylabel('|b1| (gauss)'); 
subplot(4,1,2),plot(tt,abs(gx)); ylabel('Gx (g/cm)'); 
subplot(4,1,3),plot(tt,abs(gy)); ylabel('Gy (g/cm)'); 
subplot(4,1,4),plot(tt,abs(gz)); ylabel('Gz (g/cm)'); ylim([-5,5]); 
xlabel('time (msec)'); 


end




