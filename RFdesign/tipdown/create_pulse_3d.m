function [b1,gx,gy,gz] = create_pulse_3d(x_est,n_pulse_seg,ncoils,kx,ky,kz,pulse_params,dx,dy,dz,sx,sy,sz,ishard)
% b1 is in tesla

%create RF pulses
gam=26751;                      % rad/sec/g
gambar = gam/2/pi;              % Hz/g
dt = pulse_params.dt;                 % 4usec
slew_rate = pulse_params.slew_rate;              % g/cm/sec, max slew rate
g_max = pulse_params.g_max;                      % 4g, maximum gradient amp.
%kz_area = pulse_params.kz_area;
tip_angle = pulse_params.tip_angle;
% slice_thickness = pulse_params.slice_thickness;

a_kx = kx(1:n_pulse_seg)/sx;
a_ky = ky(1:n_pulse_seg)/sy;
a_kz = kz(1:n_pulse_seg)/sz;
weights = zeros(n_pulse_seg,ncoils);
weights = x_est;
%{
for j = 1:ncoils
    weight_j = squeeze(x_est(:,:,j));
    for i = 1:n_pulse_seg
        weights(i,j) = weight_j(ky(i)+floor(dy/2)+1,kx(i)+floor(dx/2)+1);
    end
end
%}

max_pe_jump = pulse_params.max_pe_jump; %15/sx; % Interesting

weights = weights*(-1i);  %before this step, weights = b1*1i*gamma*dt, but why multiply by -i not i here? 
if ishard ~= 1
 [b1,gx,gy,gz]=full_pulse_gen_sinc(tip_angle,slice_thickness,a_kx,a_ky,max_pe_jump,kz_area,weights,ncoils,slew_rate,g_max,dt);
else
 [b1,gx,gy,gz]=full_pulse_gen_hard_3d(tip_angle,a_kx,a_ky,a_kz,max_pe_jump,weights,ncoils,slew_rate,g_max,dt);
end 

return;
