% set paths
%dosetup;

% design RF pulse and save to 'rfpulse.mat'
%main_ktCont; 

% load pulse and pad with zeros at beginning and end (TOPPE requires this)
load rfpulse;
b1 = [0; b1; 0];
gx = [0; gx; 0];
gy = [0; gy; 0];
gz = [0; gz; 0];

% set scanner hardware specs
sys = toppe.systemspecs('maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
	'maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
	'raster', 4e-6);

% write RF pulse to TOPPE file (module)
toppe.writemod('rf', b1, 'gx', gx, 'gy', gy, 'gz', gz, 'system', sys, 'ofname', 'tipdown.mod');
