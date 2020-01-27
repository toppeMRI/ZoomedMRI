function [b1,gx,gy,gz] = spokerf(d,fieldmap,mask,tipangle,n_pe,slice_thickness,ishard,xyzRange)
% code used in Jon's stfr-fmri. Copied from xsense@brooks:~/autorfdesign
% folder
%addpath spoke_daehyun;
load single_coil_tx_sensitivity;
sens = ones(size(sens));  % in the spins trajectory experiment, I used all ones as sens
%[d,sens,mask,density,dx,dy,sx,sy,ncoils,fieldmap,nfreq,pulse_params] = loadup_simulation_spoke(2, 1, 1);
%n_pe = 10; 

mask = logical(mask); 


ncand = 10;
type = 2;

ncoils = 1; !!!
lambda = 5; % penalty parameter for spoke pulse
xfov = 24;
yfov = 24;
zfov = 0.5;

dx = 64;
dy = 64;
dz = 1;

% pulse parameters for spoke pulse design
pulse_params.kz_area = 6/slice_thickness;            %-4 ~ 4 cycle/cm : 4/slice_thickness or more recommended for a good slice profile
pulse_params.tip_angle = tipangle;         %10 degree
pulse_params.slice_thickness = slice_thickness;  %
pulse_params.dt = 4*(10^-6);         %4usec
pulse_params.slew_rate = 15000;      %g/cm/sec, max slew rate
pulse_params.g_max = 4;              % 4g, maximum gradient amp.



[nrmse,kx,ky] = determine_pe_field_greedy3(d,sens,mask,n_pe,type,ncand,fieldmap,xfov,yfov,pulse_params,0,ishard);

[x_est,projection,pp_nrmse] = determine_weights_field2(xfov,yfov,kx,ky,d,mask,sens,fieldmap,pulse_params,lambda,0,ishard);

[b1,gx,gy,gz] = create_pulse(x_est,n_pe,ncoils,kx,ky,pulse_params,dx,dy,xfov,yfov,ishard);
b1 = b1*1e4; % convert to gauss
disp(sprintf('the net pulse length %f msec',size(b1,1)*pulse_params.dt*1000));

% x_range = ([0:1:dx-1]-floor(dx/2))*xfov/dx;
% y_range = ([0:1:dy-1]-floor(dy/2))*yfov/dy;
% z_range = [0];

gxtmp = gx; 
gx = gy; 
gy = gxtmp; %because daehyun use meshgrid

sens = permute(sens,[3,1,2]);  
% for i = 1:length(z_range)
%     M_pr(:,:,i)=parallel_blochsim_field(0,b1,gx,gy,gz,sens,x_range,y_range,z_range(i),pulse_params.dt,fieldmap,lambda);
% end
T1 = 1000; T2 = 100; dt = 4e-6;
M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt,fieldmap,mask,T1/1000,T2/1000);

%M_pr_tr = sum(M_pr,3);
figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default;
figure,im(abs(M_pr));colormap default;


end




