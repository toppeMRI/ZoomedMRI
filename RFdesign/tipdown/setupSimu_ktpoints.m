% set up simulations for ktpoints

%  Hao Sun : May10,2013; the university of Michigan
clear; 

simumode = 1; % 1: simulation  2: using real data do experiment
PTXtype = 1; %1:single coil; 2: parallel transmit
doSteadyState = 0; % do steady state simulation if 1;
doroi = 0;          % draw ROI (1) or use the previous one (0)
dorectangle = 0;    % tip-up only within a rectangular area in the center
plane = 'axi';  % 'sag' 'axi'
reverseK = 0;

% Tread = 3.52; %for 192 x freq encoding without fatsat puls
Tread = 4.89; %for 256 x freq encoding 62.5BW
% Tread = 5.7; %Bandwidth = 125; OPTE = 2.828ms; there is 2.12 ms before z phase encoding, what's that? !!
% Tread = 3.82; %for 192 x freq encoding without fatsat pulse 5+7.412-7.95; %ms

T1 = 520; T2 = 50; %gel phantom
T1 = 1100; T2 = 90; %gray matter
% TE = 2.2;
dt = 4e-3;             % msec

nom_fa = 16;


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
    datdir = '.';
    target_x = 1:N;
    target_y = 1:N;
    target_z = [1:5:80];
    simulated_slice = 1:length(target_z); %round(size(b0map,3)/2);
    load(sprintf('%s/Psub1.mat',datdir));
    allm = m;
    figure,im(allm);
    
else
    %     datdir = '~/lab/data/spins/sep9';
    datdir = '../data/2013-01-03';
    N = 64; %GE!!
    Nz = 65;
    FOV = 24; %cm
    FOVz = 32;
    %     N = 64; %human
    %     Nz = 49;
    %     FOV = 24; %cm
    %     FOVz = 24;
    %     im3d = loadim(sprintf('%s/P23040.7',datdir),64,32);
    %     figure,im(im3d);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P27136.7',datdir); sprintf('%s/P27648.7',datdir)],N,Nz);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P34304.7',datdir); sprintf('%s/P34816.7',datdir)],N,Nz);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P01536.7',datdir); sprintf('%s/P02048.7',datdir)],N,Nz);
    [fieldmap,allm] = measureb0([sprintf('%s/P04608.7',datdir); sprintf('%s/P05632.7',datdir)],N,Nz);
    %     [fieldmap,allm] = measureb0([sprintf('%s/P11264.7',datdir); sprintf('%s/P11776.7',datdir)],N,Nz); %GE!!
    figure,im(allm(:,:,:)); title('im');
    Nzc = round((Nz+1)/2); 
    target_z = [Nzc-2:1:Nzc+2]; %!! the targetz has to be centered to 0 because of GMRI_sense
    %     target_z = [32:1:35]; GE!!
    simulated_slice = 1:length(target_z);
    target_x = 1:N;
    target_y = 1:N;
end
b0map = fieldmap(:,:,target_z);
m = allm(:,:,target_z);
b0map = smoothim(b0map,3,0.5);
b0map = 1*b0map;

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
        if 1
            eroStr = ones(4,4);
            for iero = 1:size(b0map,3)
                %                 mask(:,:,iero) = imerode(squeeze(mask(:,:,iero)),eroStr);
                mask(:,:,iero) = imopen(squeeze(mask(:,:,iero)),eroStr);
            end
        end
        roi = mask;
        roiAll = autoROI(fieldmap,allm);
        %         load(sprintf('%s/roi',datdir));
        
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
end
sens = sens/max(abs(sens(:)));
save sensitivities sens
sens_simu = sens(:,:,:,simulated_slice);
save sensitivities_2d sens_simu
ncoils = size(sens,1);

%% define design x,y,z range for tip up
if simumode == 1
    dim = size(sens,2);
    dimz = length(target_z); 
    reso = FOV/N;
    resoz = FOVz/(Nz-1);
    xyzRangeAll.x = -FOV/2+reso/2:reso:FOV/2-reso/2;
    xyzRangeAll.y = -FOV/2+reso/2:reso:FOV/2-reso/2;
    xyzRangeAll.z = -FOVz/2:resoz:FOVz/2;%

    xyzRange.x = xyzRangeAll.x(target_x);
    xyzRange.y = xyzRangeAll.y(target_y);
    xyzRange.z = xyzRangeAll.z(target_z);
    xyzRange.z = -(xyzRange.z(end) + xyzRange.z(1))/2 + xyzRange.z; % shift the center to 0  
else
    reso = FOV/N;
    resoz = FOVz/(Nz-1);
    xyzRangeAll.x = -FOV/2+reso/2:reso:FOV/2-reso/2;
    xyzRangeAll.y = -FOV/2+reso/2:reso:FOV/2-reso/2;
    xyzRangeAll.z = -FOVz/2:resoz:FOVz/2;%
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

save designParas;
