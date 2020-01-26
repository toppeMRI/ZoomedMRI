function seq = getparams
% Define acquisition parameters

if 0
	% use Rx from SlicePlanner GUI
	roi = toppe.getroi('ROI.h5', 1); 
	seq.rotmat = roi.rotmat;
	resLoc = 0.1; % Voxel size of localizer volume (assumed to be isotropic) (cm).            
	seq.sliceOffset = resLoc*roi.scanPlaneToIsocenterDistance;      % cm
else
	% scan at iso-center, non-oblique
	seq.rotmat = eye(3);
	seq.sliceOffset = 0;   % cm
end

% Hardware limits
% NB! If passing 'sys' to writemod.m, 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' can always be a design choice, i.e., it can be at or below the physical system limit.
seq.sys = toppe.systemspecs('maxSlew', 12, 'slewUnit', 'Gauss/cm/ms', 'maxGrad', 5, 'gradUnit', 'Gauss/cm'); 

% Resolution and FOV
seq.res = [0.1 0.1 0.1];       % Resolution of high-res B0 mapping sequences (cartesian and spiral) (cm)
targetfov  = [22 22 5];        % actual FOV can be slightly different due to rounding to even matrix size
seq.matrix = [2*round(targetfov(1:2)./seq.res(1:2)/2) round(targetfov(3)/seq.res(3)];
seq.fov = seq.res .* seq.matrix;

% spiral readout parameters
seq.nLeafs = 40;              % number of spiral leafs

% fat sat pulse parameters
seq.fatsat.tbw = 1.5;
seq.fatsat.dur = 3;              % pulse duration (msec)
seq.fatsat.freqOffset = -440;    % Hz (fat frequency offset at 3T)

%
seq.rf_spoil_seed = 117;      % RF spoiling (linear phase increment) factor (degrees)
seq.fmri.nframes = 10;        % Number of time-frames.
seq.ndisdaq = 10;        % Throw away this many frames at beginning to reach steady-state

% Fat saturation pulse (bandwidth = 500 Hz)
