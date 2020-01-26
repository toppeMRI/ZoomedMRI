function main
% Create stack-of-spirals TOPPE sequence.
% This script creates 'readout.mod' and 'spoiler.mod'
% It is assumed the excitation files 'tipdown.mod' and 'tipup.mod' 
% already reside in this folder

% Usage:
%  >> addpath ..
%  >> seq = getparams();
%  >> main(seq);

%% Create modules.txt
% Entries are tab-separated.
modFileText = ['' ...
'Total number of unique cores\n' ...
'3\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'readout.mod	0	0	1\n' ...
'fatsat.mod	0	1	0\n' ...
'tipdown.mod	0	1	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);


%% Create (balanced) stack-of-spirals readout leaf. Isotropic resolution.

[g,roInfo] = toppe.utils.spiral.makesosreadout(seq.fmri.fov, seq.fmri.matrix, ...
	seq.fmri.nLeafs, 0.99*seq.sys.maxSlew, 'system', seq.sys, 'ofname', 'readout.mod', ...
	'rewDerate', 0.7, 'inout', 'in');
readoutDur = seq.sys.raster*1e3*length(roInfo.sampWin);  % msec
sampWin = roInfo.sampWin;


%% Create slab-selective excitation (tipdown.mod). Include (PRESTO) spoiler gradients.

[rf,gex,tipdown.freq] = toppe.utils.rf.makeslr(seq.fmri.flip, seq.fmri.slabThick, seq.rf.tbw, seq.rf.dur, 0, ...
                                  'ftype', seq.rf.ftype, 'system', seq.sys, 'writeModFile', false);

% Create spoiler (PRESTO) gradients.
% Will be placed on two axes for best (RF) spoiling.
nCyclesSpoil = 1;   % Gives near-optimal temporal SNR, and not so large (see one of my ISMRM abstracts, 2017 I think) 
gspoil1 = toppe.utils.makecrusher(seq.fmri.nCyclesSpoil, seq.fmri.res(3), 0, 0.7*seq.sys.maxSlew, seq.sys.maxGrad);
gspoil2 = toppe.utils.makecrusher(2*seq.fmri.nCyclesSpoil, seq.fmri.res(3), 0, 0.5*seq.sys.maxSlew, seq.sys.maxGrad);

% create tipdown.mod
rf = [0*gspoil2(:); zeros(2,1); rf(:); 0*gspoil1(:)];
gx = [1*gspoil2(:); zeros(2,1); 0*gex(:);  -gspoil1(:)];
gy = [0*gspoil2(:); zeros(2,1); 0*gex(:); 0*gspoil1(:)];
gz = [1*gspoil2(:); zeros(2,1); gex(:);  -gspoil1(:)];
rf = toppe.utils.makeGElength(rf);
gx = toppe.utils.makeGElength(gx);
gy = toppe.utils.makeGElength(gy);
gz = toppe.utils.makeGElength(gz);
toppe.writemod('rf', rf, 'gx', gx, 'gy', gy, 'gz', gz, 'ofname', 'tipdown.mod', ...
               'desc', 'RF slab excitation with PRESTO gradients', 'system', seq.sys);

%% Create fat saturation pulse
% bw = 500 Hz. Frequency offset (-440 Hz) is set in scanloop.txt.
seq.fatsat.flip = 50;
slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse
seq.fatsat.tbw = 1.5;
seq.fatsat.dur = 3;            % pulse duration (msec)
toppe.utils.rf.makeslr(seq.fatsat.flip, seq.fatsat.slThick, seq.fatsat.tbw, seq.fatsat.dur, 1e-8, ...
                       'ftype', 'ls', 'type', 'ex', 'ofname', 'fatsat.mod', 'system', seq.sys);

%% Create scanloop.txt
nz = seq.fmri.matrix(3);

% fully sampled kz sampling pattern
for ii = 1:nz
	kzFull(ii) = ((ii-1+0.5)-nz/2)/(nz/2);    % scaling is (-1 1)
end

% undersampled kz sampling pattern
if 0
	% Variable-density (non-Cartesian) kz undersampling. May be useful later.
	Rz = nz/nz_samp;               % kz acceleration factor
	kz = vardenskz(nz,Rz,3.3);    % Fully sampled center with quadratically increasing FOV outside center. Last argument is FOV(kz=0)/FOV(kzmax). 
	a_gz_max = abs((0.5-nz/2)/(nz/2));
	kzU = kz*a_gz_max;     % scaling is (-1 1)
else
	% Cartesian variable-density kz undersampling
	%load zInd;
	%kzU = kzFull(logical(zInd));
end

rfphs = 0;              % radians
rfphsLast = rfphs;
daqphs = 0;
rf_spoil_seed_cnt = 0;

% loop over (undersampled) time frames and fill in scanloop.txt,
% and write k-space values to file
if seq.writeKspace
	[rf,gx,gy] = toppe.readmod('readout.mod');  % one spiral leaf
	[kx1,ky1] = toppe.utils.g2k([gx(:,1) gy(:,1)],1);
	k1 = complex(kx1,ky1);  % single spiral 'prototype'
	ndat = size(k1,1);
	necho = 1;
	ksp.kx = NaN*ones(ndat, nz, necho, seq.fmri.nframes);   % to match the 'slice-echo-view' order of 'dat' array returned by toppe.utils.loadpfile
	ksp.ky = ksp.kx;
	ksp.kz = ksp.kx;
end

fprintf('Writing scanloop.txt for fMRI sequence\n');

toppe.write2loop('setup', 'version', 3);

for iframe = (-seq.fmri.ndisdaq+1):seq.fmri.nframes
	if ~mod(iframe,10)
		fprintf([repmat('\b',1,20) sprintf('%d of %d', iframe+seq.fmri.ndisdaq, seq.fmri.nframes+seq.fmri.ndisdaq )]);
	end

	% Set kz sampling pattern for this frame.
	kz1 = kzFull;

	% set 'view' data storage index
	if iframe < 1 
		dabmode = 'off';
	else
		dabmode = 'on';
	end

	for iz = 1:nz

		a_gz = ((iz-1+0.5)-nz/2)/(nz/2);   % z phase-encode amplitude, scaled to (-1,1) range

		for ileaf = 1:seq.fmri.nLeafs
			if seq.fmri.doFatsat
  	 			toppe.write2loop('fatsat.mod', 'RFphase', rfphs, 'Gamplitude', [0 0 0]', ...
					'RFoffset', seq.fatsat.freqOffset);
				% rfphs = rfphs + (seq.fmri.rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
				% rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
			end

   		% rf excitation module (includes PRESTO gradients)
   		toppe.write2loop('tipdown.mod', 'RFphase', rfphs, ...
					'rotmat',   seq.rotmat, ...
					'RFoffset', round(tipdown.freq));

   		% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
			slice = iz;
			echo = 1;
			view = (max(iframe,1)-1)*seq.fmri.nLeafs + ileaf;
			phi = 2*pi*(ileaf-1)/seq.fmri.nLeafs;                 % leaf rotation angle (radians)
   		toppe.write2loop('readout.mod', 'DAQphase', rfphsLast, ...
				'slice', slice, 'echo', echo, 'view', view, ...
				'rotmat', seq.rotmat, ...        % prescribed scan plane rotation (logical frame)
				'rot',    phi, ...               % in-plane rotation (in logical frame)
				'dabmode', dabmode, 'Gamplitude', [1 1 a_gz]');

	   	% update rf phase (RF spoiling)
			rfphsLast = rfphs;
			rfphs = rfphs + (seq.fmri.rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;

			% kspace info for this TR
			if strcmp(dabmode, 'on') & seq.writeKspace
				k1tmp = k1.*exp(1i*phi)*seq.fmri.fov(1)/seq.fmri.matrix(1);         % convert from cycles/cm to cycles/sample
				ksp.kx(:,slice,echo,view) = real(k1tmp);
				ksp.ky(:,slice,echo,view) = imag(k1tmp);
				ksp.kz(:,slice,echo,view) = kz1(iz)/2;         % cycles/sample
			end

			ksp.kinfo(slice,echo,view).ileaf = ileaf;     % might come in handy
			ksp.kinfo(slice,echo,view).rot = phi;         % this too
		end
	end
end
fprintf('\n');
toppe.write2loop('finish');

%% Save k-space trajectory and other sequence info 
if seq.writeKspace
	ksp.kx = single(ksp.kx);
	ksp.ky = single(ksp.ky);
	ksp.kz = single(ksp.kz);
	fprintf('Writing ksp.mat...');
	%save -v7.3 ksp ksp
	save ksp ksp
	fprintf('done\n');
	kspfile = 'ksp.mat';
else
	kspfile = '';
end

%% create tar file
system(sprintf("tar czf scan,fmri.tgz ../getparams.m *.m *.mod modules.txt scanloop.txt %s", kspfile));

fprintf('Scan time for 3D spiral fmri sequence: %.2f min\n', toppe.getscantime/60);

%% display sequence
%toppe.playseq(2, 'tpause', 0.1);

%% convert to Pulseq
if 0
% addpath ~/gitlab/toppe/pulseq/
cd tar
ge2seq('scan.tar');
seq = mr.Sequence();
seq.read('out.seq');
seq.plot_sg('TimeRange', [10 10.05]);
cd ..
end


