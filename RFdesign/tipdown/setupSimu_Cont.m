% set up simulations for ktpoints

%  Hao Sun : May10,2013; the university of Michigan

% field map location
fmapdir = '~/gitlab/zoomed-mri/b0map/';
save fmappath fmapdir

% these are set in /usr/g/bin/recon812 on scanner
%exptype = 'slr';   % 'iv': use tailored pulse designed here; 'slr': use SLR tipdown pulse in ~/projects/myelin/tipdown.mod
exptype = 'iv';
doroi = false;           % draw mask (1) or use the previous one (0)

if strcmp(exptype,'slr')
	Tread = 7.56;   % ms. For STFR myelin water fraction mapping
else
	Tread = 16; % ms. For IV STFR fMRI
end
save Tread Tread

simumode = 2; % 1: simulation  2: using real data do experiment
PTXtype = 1; %1:single coil; 2: parallel transmit
doSteadyState = 0; % do steady state simulation if 1;
dorectangle = 0;    % tip-up only within a rectangular area in the center
plane = 'axi';  % 'sag' 'axi'
reverseK = 0;

% set these parameters to match your experiment
T1 = 520; T2 = 50; %gel phantom
T1 = 1300; T2 = 50; % T2star of grey matter
%T1 = inf; T2 = inf;
% TE = 2.2;

dt = 4e-3;             % msec

zthick = 0.125; % cm
ncycles = 1.5; % number of gradient crusher cycles in a zthick

pseqtype = 0;
if pseqtype==1               % bSSFP
    isbssfp = 1;
else
    isbssfp = 0;             % STFR
end
if pseqtype==2
    tipdownonly = 1;
else
    tipdownonly = 0;
end


%%
%% load fieldmap and m (magnitude image)
%%
datdir = '.'; 
N = 64; %GE!!
Nz = 24;
FOV = 24; %cm
FOVz =24;
  
load(sprintf('%s/fieldmap.mat',fmapdir));
%fieldmap = -fieldmap;   % for testing
allm = m;
figure,im(allm(:,:,:)); title('im');
%target_z = 7:2:Nz-2; % for ktraj_paper
target_z = 1:1:Nz; % used for reduced FOV design, since you want to cover the whole object
simulated_slice = 1:length(target_z);
target_x = 1:N;
target_y = 1:N;
save targets target_x target_y target_z

b0map = fieldmap(:,:,target_z);
m = allm(:,:,target_z);
b0map = smoothim(b0map,3,0.5);
fieldmap = smoothim(fieldmap,3,0.5);
%b0map = 1*b0map; !!!

b0map_simu = b0map(:,:,simulated_slice);
save b0map b0map
save fieldmap fieldmap
save b0map_simu b0map_simu
save curim m;


%%
%% define inner-volume (IV.mat)
%%
if doroi  % create IV.mat
	%figure,im(m);
	roipicker(m);     % saves ROI in IV.mat
	h = figure(3);    % move figure 3 to front. Figure 3 is the axial view in roipicker.
	waitfor(h);       % continue execution after Figure 3 is closed. 
	%title('Draw ROI');
else 
	% load existing IV.mat
end
%load ../../ivRx/IV;
load IV

%%
%% define mask (care-region)
%%  For legacy code reasons, roi, roifull, and maskfull are also defined (= mask)
%%
mask = autoROI(b0map,m);
maskfull = autoROI(fieldmap,allm);
if 1 % some noise pixels can increase the target BW a lot
	eroStr = ones(4,4);
	for iero = 1:size(b0map,3)
		% mask(:,:,iero) = imerode(squeeze(mask(:,:,iero)),eroStr);
		mask(:,:,iero) = imopen(squeeze(mask(:,:,iero)),eroStr);
	end
	for iero = 1:size(fieldmap,3)
		maskfull(:,:,iero) = imopen(squeeze(maskfull(:,:,iero)),eroStr);
	end
end
roifull = maskfull;
save roifull roifull;
roi = mask; 
save mask mask;
save roi roi;
roi_simu = roi(:,:,simulated_slice);
save roi_simu roi_simu;
figure,im(roi);
figure,im(fieldmap(:,:,target_z).*roi,'cbar'); colormap('default');


%%
%% load coil sensitivities, and create sensitivities.mat
%%
if PTXtype ==2
    load sensitivities_8coil.mat;
    for c = 1:size(sens,3)
        dimz = length(target_z);
        tmp = repmat(sens(:,:,c),[1,1,dimz]);
        sensall(c,:,:,:) = tmp;
    end
    sens = sensall;
else
    sens(1,:,:,:) = ones(size(m));
    if simumode == 2 % becareful when use this sens map, since some region is 0.
        %       sens(1,:,:,:) = abs(allm(:,:,target_z))/max(max(max(abs(allm(:,:,target_z)))));
    end
end
sens = sens/max(abs(sens(:)));
save sensitivities sens
sens_simu = sens(:,:,:,simulated_slice);
save sensitivities_2d sens_simu
ncoils = size(sens,1);


%%
%% define design x,y,z range
%%
reso = FOV/N;
resoz = FOVz/Nz; %(Nz-1); % This should be right according to (1) stfr5.e (2) agree with simu
xyzRangeAll.x = -FOV/2:reso:(FOV/2-reso);
xyzRangeAll.y = -FOV/2:reso:(FOV/2-reso);
xyzRangeAll.z = -FOVz/2:resoz:(FOVz/2-resoz); % this agree with the scanner acquisition
xyzRange.x = xyzRangeAll.x(target_x);
xyzRange.y = xyzRangeAll.y(target_y);
xyzRange.z = xyzRangeAll.z(target_z);

simuRange.x = xyzRange.x;
simuRange.y = xyzRange.y;
simuRange.z = xyzRange.z(simulated_slice);

save Range xyzRange xyzRangeAll simuRange

save designParas;
