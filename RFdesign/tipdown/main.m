% set paths
dosetup;

% design RF pulse and save to 'rfpulse.mat'
%main_ktCont; 

% load pulse
load rfpulse;

% set scanner hardware specs
sys = toppe.systemspecs('maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
	'maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
	'raster', 4e-6);

% try to write to file
% toppe.writemod('rf', b1, 'gx', gx, 'gy', gy, 'gz', gz, 'system', sys, 'ofname', 'tipdown.mod');
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
toppe.writemod('rf', b1, 'gx', gx, 'gy', gy, 'gz', gz, 'system', sys, 'ofname', 'tipdown.mod');

% display the .mod file
% Does it satisfy PNS limits?
toppe.plotmod('tipdown.mod');
