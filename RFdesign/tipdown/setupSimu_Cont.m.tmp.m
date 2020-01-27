% set up simulations for ktpoints

%  Hao Sun : May10,2013; the university of Michigan
% clear;

simumode = 2; % 1: simulation  2: using real data do experiment
pregen  = 1;  % for pre-generated rFOV pulse: use a large cylinder. Fieldmap is set to zero.
PTXtype = 1; %1:single coil; 2: parallel transmit
doSteadyState = 0; % do steady state simulation if 1;
doroi = 0;          % draw ROI (1) or use the previous one (0)
dorectangle = 0;    % tip-up only within a rectangular area in the center
plane = 'axi';  % 'sag' 'axi'
reverseK = 0;

% Tread = 3.52; %for 192 x freq encoding without fatsat puls
Tread = 4.89/2; %for 256 x freq encoding 62.5BW
% Tread = 5.7; %Bandwidth = 125; OPTE = 2.828ms; there is 2.12 ms before z phase encoding, what's that? !!
% Tread = 3.82; %for 192 x freq encoding without fatsat pulse 5+7.412-7.95; %ms


T1 = 520; T2 = 50; %gel phantom
T1 = 1100; T2 = 90; %gray matter
% TE = 2.2;
dt = 4e-3;             % msec


zthick = 0.125; % cm
ncycles = 2; % number of gradient crusher cycles in a zthick

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





%% load B0 map
% [b0map,m] = measureb02d([sprintf('%s/Pte1.7',datdir); sprintf('%s/Pte2.7',datdir)]);     % Hz, size 64x64
% [b0map,m] = measureb0([sprintf('%s/Pte13d.7',datdir); sprintf('%s/Pte23d.7',datdir)]);     % Hz, size 64x64
if simumode == 1
    FOV = 24; % cm
    FOVz = 8; % cm
    Nz = 80;
    N = 64;
    datdir = '..';
    target_x = 1:N;
    target_y = 1:N;
    target_z = [40:5:75];
    simulated_slice = 1:length(target_z); %round(size(b0map,3)/2);
    load(sprintf('%s/Psub1.mat',datdir));
    allm = m;
    figure,im(allm);
    
else
    
    %     datdir = '~/lab/data/spins/sep9';
    %     datdir = '../data/2013-01-03';
    %     datdir = '~/lab/data/trajDesign/nov1/';
    datdir = '.'; 
    N = 64; %GE!!
    Nz = 24;
    FOV = 24; %cm
    FOVz =24;
    %     N = 64; %human
    %     Nz = 49;
    %     FOV = 24; %cm
    %     FOVz = 24;
    %     im3d = loadim(sprintf('%s/P23040.7',datdir),64,32);
    %     figure,im(im3d);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P27136.7',datdir); sprintf('%s/P27648.7',datdir)],N,Nz);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P34304.7',datdir); sprintf('%s/P34816.7',datdir)],N,Nz);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P01536.7',datdir); sprintf('%s/P02048.7',datdir)],N,Nz);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P53760.7',datdir); sprintf('%s/P54784.7',datdir)],N,Nz);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P08192.7',datdir);
    %     sprintf('%s/P08704.7',datdir)],N,Nz);% oct 22
    %     [fieldmap,allm] = measureb0([sprintf('%s/P15360.7',datdir); sprintf('%s/P15872.7',datdir)],N,Nz);
    
    % the following is a gel phantom b0map used for simu/testing
    [fieldmap,allm] = measureb0(['P07680.7'; 'P09216.7'],N,Nz);
        
%     [fieldmap,allm] = measureb0_sos3d('readout.wav', ['P49152.7'; 'P50176.7'],1);

%    load ims_fmap; 
%    OPTE = [3, 5.3] * 1e-3; %sec
%    fieldmap = angle(ims2./ims1)/(OPTE(2)-OPTE(1))/(2*pi);    % Hz
%    allm = ims1./max(abs(ims1(:)));

%        datdir = '~/lab/data/rFOV/oct22/'; %for Oct23 scan jon
%        [fieldmap, allm] = measureb0_sos3d([datdir 'readoutOct22b0map.wav'], [[datdir,'P63488.7']; [datdir,'P64000.7']],1);
%        Nz = size(fieldmap,3);
%        FOVz = 18;
%        for iz = 1:Nz
%        fieldmapSmall(:,:,iz) = imresize(fieldmap(:,:,iz), [N,N]);
%        allmSmall(:,:,iz) = imresize(allm(:,:,iz), [N,N]);
%        end
%        fieldmap = fieldmapSmall;
%        allm = allmSmall;

%    datdir = '~/lab/data/trajDesign/nov7_fieldmap_fromSydney/';
%    [fieldmap, allm] = measureb0_sos3d([datdir 'readout.wav'], [[datdir,'P13312.7']; [datdir, 'P13824.7']],1); 
%    Nz = size(fieldmap,3); 
%    FOVz = 18; 

    if 0
        fieldmap = zeros(size(fieldmap));
        allm = repmat(allm(:,:,round(end/2)),[1,1,size(allm,3)]);
    end
    figure,im(allm(:,:,:)); title('im');
    %     Nzc = round((Nz+1)/2);
    %     target_z = [Nzc-2:1:Nzc+2]; %!! the targetz has to be centered to 0 because of GMRI_sense
    %     target_z = [32:1:35]; GE!!
	if pregen
	    target_z = 1:Nz; % used for reduced FOV disign, since you want to cover the whole object
	else
	    target_z = 7:2:Nz-2; % for ktraj_paper
	end
    simulated_slice = 1:length(target_z);
    target_x = 1:N;
    target_y = 1:N;
end
b0map = fieldmap(:,:,target_z);
m = allm(:,:,target_z);
b0map = smoothim(b0map,3,0.5);
b0map = 1*b0map; !!!

b0map_simu = b0map(:,:,simulated_slice);
save b0map b0map
save b0map_simu b0map_simu
save curim m;

%% define ROI
if simumode == 1
    mask = mask(:,:,target_z);
    roi = mask;
else
    if doroi
        figure,im(m);
        title('Draw ROI');
        roi = roipoly;
        save(sprintf('%s/roi',datdir),'roi');
    else
        mask = autoROI(b0map,m);
        if 1 % some noise pixels can increase the target BW a lot
            eroStr = ones(4,4);
            for iero = 1:size(b0map,3)
                %                 mask(:,:,iero) = imerode(squeeze(mask(:,:,iero)),eroStr);
                mask(:,:,iero) = imopen(squeeze(mask(:,:,iero)),eroStr);
            end
        end
        roi = mask; 
        if pregen % for pre-generated rFOV pulse: use a large cylinder. Fieldmap is set to zero.
           
           %         roi = mask;
           %         roiAll = autoROI(fieldmap,allm);
           %         load(sprintf('%s/roi',datdir));
        [Nx, Ny, Nz] = size(b0map); 
        xr = 1:Nx; yr = 1:Ny; 
        [xxtmp,yytmp] = ndgrid(xr-32,yr-32);
        radiusX = 11/24*Nx;
        radiusY = 11/24*Ny; 
        circInd = (xxtmp/radiusX).^2 + (yytmp/radiusY).^2 < 1;
        z1 = 1;
        z2 = Nz;
        circInd = repmat(circInd, [1 1 z2-z1+1]);
        roi = false(size(b0map)); 
        roi(:,:,z1:z2) = circInd;
        mask = roi; 
        b0map = zeros(size(b0map)); 
        fieldmap = zeros(size(fieldmap)); 
        end 
    end
end
save mask mask;
save roi roi;
roi_simu = roi(:,:,simulated_slice);
save roi_simu roi_simu;
figure,im(roi);
figure,im(fieldmap(:,:,target_z).*roi,'cbar'); colormap('default');

%% load coil sensitivities, and create sensitivities.mat
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

%% define design x,y,z range for tip up
if simumode == 1
    %     dim = size(sens,2);
    %     dimz = length(target_z);
    %     reso = FOV/N;
    %     resoz = FOVz/(Nz-1);
    %     xyzRangeAll.x = -FOV/2+reso/2:reso:FOV/2-reso/2;
    %     xyzRangeAll.y = -FOV/2+reso/2:reso:FOV/2-reso/2;
    %     xyzRangeAll.z = -FOVz/2:resoz:FOVz/2;%
    %
    %     xyzRange.x = xyzRangeAll.x(target_x);
    %     xyzRange.y = xyzRangeAll.y(target_y);
    %     xyzRange.z = xyzRangeAll.z(target_z);
    %     xyzRange.z = -(xyzRange.z(end) + xyzRange.z(1))/2 + xyzRange.z; % shift the center to 0
    %
    dim = size(sens,2);
    dimz = length(target_z);
    FOVz = (target_z(end) - target_z(1)+target_z(2)-target_z(1))*0.1;
    reso = FOV/dim;
    resoz = FOVz/(dimz-1); % this agree with stfr5.e acquisition, the Gmri_sense needs to be scaled for this.
    xyzRange.x = -FOV/2:FOV/dim:FOV/2-reso;
    xyzRange.y = -FOV/2:FOV/dim:FOV/2-reso;
    xyzRange.z = -FOVz/2:resoz:FOVz/2; % this agree with the scanner acquisition, the Gmri_sense needs to be shifted and scaled
else
    reso = FOV/N;
    resoz = FOVz/(Nz-1); % This should be right according to (1) stfr5.e (2) agree with simu
    xyzRangeAll.x = -FOV/2:reso:FOV/2-reso;
    xyzRangeAll.y = -FOV/2:reso:FOV/2-reso;
    xyzRangeAll.z = -FOVz/2:resoz:FOVz/2; % this agree with the scanner acquisition
    %    xyzRangeAll.z = (-Nz/2:Nz/2-1)*resoz;
    %     xyzRangeAll.z = xyzRangeAll.z - resoz/2;
    %     xyzRangeAll.z = -FOVz/2:resoz:FOVz/2-resoz;
    %     xyzRangeAll.z =  xyzRangeAll.z*1.667;
    xyzRange.x = xyzRangeAll.x(target_x);
    xyzRange.y = xyzRangeAll.y(target_y);
    xyzRange.z = xyzRangeAll.z(target_z);
end
simuRange.x = xyzRange.x;
simuRange.y = xyzRange.y;
simuRange.z = xyzRange.z(simulated_slice);

%% define weighting matrix
% [gbx gby gbz] = gradient(b0map);
% gb = gbx.^2+gby.^2 +gbz.^2;
gb = ones(size(b0map));
%gb = smoothim(gb, 5, 2);
%gb(:,:,1:6) = 0;
%gb(:,:,11:end) = 0;
weightIm = (gb).^0.5; % same dimension as image
% weight = sparse(diag(weightIm(roi)));
weight = diag_sp(weightIm(roi));

save weighting weight weightIm gb;

save designParas;
