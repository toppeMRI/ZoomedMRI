function main
% Create dynamic 3D stack-of-spirals, RF-spoiled STFR sequence in TOPPE format.
%
% This script creates the following files:
%  'modules.txt'
%  'readout.mod' 
%  'fatsat.mod' 
%  'spoiler.mod' 
%  'scanloop.txt'
%
% It is assumed that the excitation module files 'tipdown.mod' and 'tipup.mod' 
% have already been created and reside in this folder.
% You can design those waveforms any way you like, then use 'toppe.writemod'
% to write them to .mod files.

seq = getparams;

%% Create modules.txt
% File lists the various .mod files that are to be referenced in scanloop.otxt
% Entries are tab-separated.
modFileText = ['' ...
'Total number of unique cores\n' ...
'5\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'fatsat.mod	0	1	0\n' ...
'spoiler.mod	0	0	0\n' ...
'tipdown.mod	0	1	0\n' ...
'readout.mod	0	0	1\n' ...
'tipup.mod	0	1	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);


%% Create readout.mod. (balanced) stack-of-spirals readout leaf. Isotropic resolution.
[g,roInfo] = toppe.utils.spiral.makesosreadout(seq.fov, seq.matrix, seq.nLeafs, seq.sys.maxSlew, ...
	'ofname', 'readout.mod', ...
	'system', seq.sys, ...
	'inout',  seq.spiralType);
%	'rewDerate', 0.7, ...
readoutDur = seq.sys.raster*1e3*length(roInfo.sampWin);  % msec
sampWin = roInfo.sampWin;


%% Create fatsat.mod 
slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse
flip = 90;
toppe.utils.rf.makeslr(flip, slThick, seq.fatsat.tbw, seq.fatsat.dur, 1e-8, ...
                       'ftype', 'ls', 'type', 'ex', 'ofname', 'fatsat.mod', 'system', seq.sys);

%% Create spoiler.mod
%gspoil = toppe.utils.makecrusher(seq.nCyclesSpoil, seq.res(1), 0, 0.7*seq.sys.maxSlew/sqrt(2), seq.sys.maxGrad/sqrt(2));
gspoil = toppe.utils.makecrusher(seq.nCyclesSpoil, seq.res(1), 0, seq.sys.maxSlew, seq.sys.maxGrad);
toppe.writemod('ofname', 'spoiler.mod', 'gx', gspoil, 'gz', gspoil, 'system', seq.sys);


%% Create scanloop.txt

% define various loop variables
nz = seq.matrix(3);        % number of kz-encodes
rfphs = 0;                 % RF phase variable (to be updated in loop) (radians)
rfSpoilSeedCnt = 0;        % RF 'shot' counter

% scan loop (dynamic 3D stack-of-spirals)
toppe.write2loop('setup', 'version', 3);                    % Initialize 'scanloop.txt'
for iframe = 1:seq.nframes                                  % temporal frame index
	for iz = 1:nz
		a_gz = ((iz-1+0.5)-nz/2)/(nz/2);                      % z phase-encode amplitude, scaled to (-1,1) range

		for ileaf = 1:seq.nLeafs
			% Fat saturation
  	 		toppe.write2loop('fatsat.mod', 'RFphase', rfphs, ...
				'Gamplitude', [0 0 0]', ...                        % Turn off any gradients that may exist in fatsat.mod
				'RFoffset', seq.fatsat.freqOffset);                % -440 Hz at 3T
  	 		toppe.write2loop('spoiler.mod');                      % Crush the fat we just excited

   		% Inner-volume excitation
   		toppe.write2loop('tipdown.mod', 'RFphase', rfphs);

   		% Data acquisition (readout). 
			% Data is stored in 'slice', 'echo', and 'view' indeces (GE-specific).
			% We 'cram' both the temporal and spiral leaf indeces into the 'view' index
			view = (max(iframe,1)-1)*seq.nLeafs + ileaf;               % A bit awkward but it works
			phi = 2*pi*(ileaf-1)/seq.nLeafs;                           % leaf rotation angle (radians)
   		toppe.write2loop('readout.mod', 'DAQphase', rfphs, ...     % Receive phase must match RF phase
				'slice', iz, 'echo', 1, 'view', view, ...
				'rotmat', seq.rotmat, ...                               % prescribed scan plane rotation (logical frame)
				'rot',    phi, ...                                      % in-plane rotation (in logical frame)
				'dabmode', 'on',...
				'Gamplitude', [1 1 a_gz]');

			% Tip-up pulse and spoiler
   		toppe.write2loop('tipup.mod', 'RFphase', rfphs);           % Tip-up phase must match tip-down phase
  	 		toppe.write2loop('spoiler.mod');                           % Crush any transverse signal left over from tip-up pulse

	   	% update rf phase (RF spoiling)
			rfphs = rfphs + (seq.rfSpoilSeed/180 * pi)*rfSpoilSeedCnt ;  % radians
			rfSpoilSeedCnt = rfSpoilSeedCnt + 1;
		end
	end
end
toppe.write2loop('finish');

%% Pre-view sequence
%toppe.playseq(2, 'tpause', 0.1);


%% Execute sequence on scanner

% GE
% Copy the .mod files, modules.txt, and scanloop.txt to /usr/g/bin/ on scanner and scan with toppev3

% Siemens
% Convert to .seq file using 'ge2seq.m' (WIP)
% Example:
% addpath ~/gitlab/toppe/pulseq/
% ge2seq('scan.tar');    % scan.tar contains the .mod files, modules.txt, and scanloop.txt
% seq = mr.Sequence();
% seq.read('out.seq');
% seq.plot_sg('TimeRange', [10 10.05]);


