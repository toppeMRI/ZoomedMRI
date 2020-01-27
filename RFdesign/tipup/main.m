
% First design a 'pre-phasing' pulse so we can
% understand what the pulse is doing.
Trf = 2e-3;               % RF pulse duration (sec)
TE  = 3e-3;               % TE (sec) (determines target excitation)
signOfTargetPhase = -1;   % 'pre-phasing'
lambda = 0.6;             % regularization paramter
wn = [-20:0.1:20];        % target frequency range
b1 = spectralRF(Trf, TE, wn, signOfTargetPhase, lambda, 'tipdown');

% Design a tipup pulse for an STFR sequence with free precession time 5e-3 s: 
Tfree = 5e-3;  % s
signOfTargetPhase = +1;   % 'post-phasing'
figure;
b1 = spectralRF(Trf, Tfree, wn, signOfTargetPhase, lambda, 'tipup');

return;

% set scanner hardware specs
sys = toppe.systemspecs('maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
	'maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
	'raster', 4e-6);

% try to write to file
% toppe.writemod('rf', b1, 'gx', gx, 'gy', gy, 'gz', gz, 'system', sys, 'ofname', 'tipup.mod');
% Fails! Because TOPPE requires all waveforms to:
% (1) start and end at zero, and (2) be on a 4-sample (16us) boundary

b1 = [0; b1; 0];
gx = [0; gx; 0];
gy = [0; gy; 0];
gz = [0; gz; 0];

nextra = 4 - mod(length(b1),4);
b1 = [b1; zeros(nextra,1)];
gx = [gx; zeros(nextra,1)];
gy = [gy; zeros(nextra,1)];
gz = [gz; zeros(nextra,1)];

% try again
toppe.writemod('rf', b1, 'gx', gx, 'gy', gy, 'gz', gz, 'system', sys, 'ofname', 'tipup.mod');

% display the .mod file
toppe.plotmod('tipup.mod');
