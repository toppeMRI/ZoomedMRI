function seq = getparams
% Define acquisition parameters for all scans (each scan resides in its own subfolder)
% 
% This file is https://gitlab.com/fMRI/toppe-sequences/stack-of-spirals-presto-bold-fmri/blob/master/getparams.m

isTest = false;           % True: Create a sequence containing only a few time frames (for testing recon, etc)
seq.writeKspace = false;   % Write k-space locations for entire scan to a (large) .mat file. Units: cycles/cm

% Scan plane prescription (rotation and slice offset)
if 0
	% use Rx from TOPPE GUI
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
seq.sys = toppe.systemspecs('maxSlew', 11, 'slewUnit', 'Gauss/cm/ms', 'maxGrad', 5, 'gradUnit', 'Gauss/cm'); 

% Resolution and FOV
seq.res.hi = [0.1 0.1 0.1];    % Resolution of high-res B0 mapping sequences (cartesian and spiral) (cm)
seq.res.lo = [0.3 0.3 0.15];   % Resolution of BOLD fMRI sequence (cm). Note high res in z (z-shimming).
targetfov  = [22 22 5];        % actual FOV can be slightly different due to rounding to even matrix size
seq.matrix.hi = [2*round(targetfov./seq.res.hi/2)];
seq.matrix.lo = [2*round(targetfov./seq.res.lo/2)];
seq.fov.hi = seq.res.hi .* seq.matrix.hi;
seq.fov.lo = seq.res.lo .* seq.matrix.lo;

% Hi-res B0 mapping sequences (Cartesian and spiral)
seq.b0.res = seq.res.hi;
seq.b0.fov = seq.fov.hi;
seq.b0.matrix = seq.matrix.hi;
seq.b0.oprbw = 125/2;             % readout bandwidth
seq.b0.nCyclesSpoil = 1;
seq.b0.flip = 4;                 % excitation angle (degrees)
seq.b0.rf_spoil_seed = 117;      % a common choice
seq.b0.TE = [0 1 2.3];           % msec
seq.b0.spiral.nLeafs = 40;
seq.b0.doFatsat = false;          % for spiral B0 mapping sequence
seq.b0.slabThick = 0.8*seq.fov.hi(3);

% Stack-of-spirals PRESTO fMRI sequence (low res)
seq.fmri.res = seq.res.lo;
seq.fmri.fov = seq.fov.lo;
seq.fmri.matrix = seq.matrix.lo;
seq.fmri.nLeafs = 3;             % Number of spiral rotations (leafs) for full k-space sampling.
seq.fmri.nCyclesSpoil = 1;       % Gives near-optimal temporal SNR, and not so large (see one of my ISMRM abstracts, 2017 I think) 
seq.fmri.flip = 10;              % excitation angle (degrees)
seq.fmri.slabThick = 0.8*seq.fov.hi(3);  % make it a bit thinner than b0 map to make sure we have good B0 map near slab edge
seq.fmri.rf_spoil_seed = 150;    % For 6/18/30/42 kz platters per frame and rf_spoil_seed=150, ...
                                 %   we have mod(nz_samp*rf_spoil_seed,360)=180 which is what we need for improved RF spoiling...
                                 %   using our method in Magn Reson Med. 2016 Jun;75(6):2388-93. doi: 10.1002/mrm.25843.
seq.fmri.nframes = 10;           % Number of time-frames.
seq.fmri.ndisdaq = 10;           % Throw away this many frames at beginning to reach steady-state
seq.fmri.doFatsat = true;

% Slab-selective excitation parameters (common to all sequences)
seq.rf.tbw = 8;                             % time-bandwidth product of SLR pulse 
seq.rf.dur = 2;                             % RF pulse duration (msec)
seq.rf.ftype = 'ls';                        % least-squares SLR pulse (another good option for 3D imaging is 'min')

% Fat saturation pulse (bandwidth = 500 Hz)
seq.fatsat.flip = 70;
seq.fatsat.slThick = 1000;       % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse
seq.fatsat.tbw = 1.5;
seq.fatsat.dur = 3;              % pulse duration (msec)
seq.fatsat.freqOffset = -440;    % Hz

