function [d,sens,mask,density,dx,dy,sx,sy,ncoils,b0map,nfreq,pulse_params] = loadup_simulation(type,slice_num,sub_num)

% load transmission sensitivity pattern(B1 field map) into sens
% size : [dy dx ncoils]
load single_coil_tx_sensitivity.mat
%load sensitivities;
% sens = cat(3, sens(:,:,1), sens(:,:,2)); % Hao

%sens = ones(64,64);

% number of coils
ncoils = size(sens,3);

% dimension in x and y (space)
% 64x64 sampling grid for x and y direction
dx = 64;
dy = 64;

% size of the field of view(FOV) in cm
sx = 24; %24 cm for x direction
sy = 24; %24 cm for y direction

%fieldmap setting
subject = sprintf('Psub%d.mat',sub_num);
load(subject);
b0map = squeeze(fieldmap(:,:,slice_num));
% load b0map_linearShim; 
%b0map = 0*b0map; % Hao: set to 0 if ignore b0 effect
% region of the interest
%[xpos, ypos] = meshgrid([0:1:dx-1],[0:1:dy-1]);
%mask = (((xpos-32).^2)+ ((ypos-32).^2))/27^2 < 1;  % a circle at the center of FOV with a radius of 27 pixel(diamter = 24/64*27 cm)
mask = squeeze(mask(:,:,slice_num));

% density
density = squeeze(m(:,:,slice_num));

% number of phase encoding locations
nfreq = 3;

%RF pulse and gradient parameter structure, pulse_params,
%contains following fields.
%kz_area : the length of a kz-line in the Echo-Volumar trajectory
%          unit : cycle/cm
%tip_angle : tip angle of the RF pulse
%            unit : degree
%slice_thickness : the thickness of the excited slice
%                  unit : cm
%dt : sampling interval for the pulse and the gradient waveforms
%     unit : sec
%slew_rate : maximum slew rate of the gradients
%     unit : g/cm/sec
%g_max : maximum amplitude of the gradients
%     unit : gauss

pulse_params.kz_area = 8;            %-4 ~ 4 cycle/cm : 4/slice_thickness or more recommended for a good slice profile
pulse_params.tip_angle = 10;         %10 degree
pulse_params.slice_thickness = 0.5;  %0.5cm
pulse_params.dt = 4*(10^-6);         %4usec
pulse_params.slew_rate = 15000;      %g/cm/sec, max slew rate
pulse_params.g_max = 4;              % 4g, maximum gradient amp.

if type == 1
    % a homogeneous excitation pattern
    d = ones(dy,dx);
elseif type == 2
    % pre-phased pattern for RF-SSFP
    d = ones(dx,dy).*exp(1i*2*pi*b0map*10*10^(-3)/2);
elseif type == 3
    d = zeros(dy,dx); 
    d(dy/2-10:dy/2+10, dx/2-10:dx/2+10) = 1; 
end


% clear the regions which are not included in the region of the interest
d(~mask)=0;
b0map(~mask) = 0;
density(~mask)=0;

for i = 1:ncoils
    sens_i = sens(:,:,i);
    sens_i(~mask)=0;
    sens(:,:,i)=sens_i;
end


return;
