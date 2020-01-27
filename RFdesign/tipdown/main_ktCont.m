%function main_ktCont(n_pe,xoff,yoff,reply)
% main file for continuous trajectory based on kt points

n_pe = 80;      % determines pulse duration
niter = 20;

xoff = 0;
yoff = 0;
reply = 'y';

extype = 3; % 1: uniform excitation; 2: pre-phasing; 3: cube/hockey puck or manual 3D IV; 4: slice selective pre-phasing;
trajtype = 0;  %0: extended KT-points; 1: SPINS; 2: SoS
dorect = 1;

%% setup
close all;
setupSimu_Cont;  !!! check the b0map is scaled or not

switch exptype
	case 'iv'
		tipangle = 17;  % for maximum grey matter signal
	case 'slr'
		tipangle = 15;  % for ~/projects/myelin/psd/
end
save tipangle tipangle

%dospin = 0 % compare with spins trajectory
%doSoS = 0; % compare with stack of spiral trajectory
isBssfp = 1; % balanced gradient
ishard = 1; % used to determine the subpulse type (spoke only) and its duration (only for the b0 correction)
tol = 0.5;

dosmooth = 0; % smooth the edge of target excitation pattern

if extype == 1
    d = ones(size(b0map_simu));
    d(~mask) = 0;
end

if extype == 2 % pre-phaseing;
    %n_pe = 70;
    if 0
        mask = logical(d);
        roi = logical(d);
        save mask mask;
        save roi roi;
        d = exp(1i*b0map_simu*Tread/1000*2*pi);
        %     d(:,:,1:ceil(size(b0map_simu,3)/2)-3) = 0;
        %     d(:,:,ceil(size(b0map_simu,3)/2)+3:end) = 0;
        d(~mask) = 0;
    else
        b0map_simu = b0map_simu; % !!! increase the field map for pre-phasing problem
        b0map = b0map;
        d = exp(1i*b0map_simu*Tread/1000*2*pi);
        d(~mask) = 0;
    end
end


%n_pe = 70; % for ktraj paper
%n_pe = 80;
%load ../../ivRx/IV
load IV
d = 1.0*IV;

if extype == 4 %slice selective pre-phasing;
    n_pe = 60;
    dp = exp(1i*b0map_simu*Tread/1000*2*pi);
    d = zeros(size(dp));
    d(:,:,ceil(size(b0map_simu,3)/2-1):ceil(size(b0map_simu,3)/2)+3) = dp(:,:,ceil(size(b0map_simu,3)/2-1):ceil(size(b0map_simu,3)/2)+3);
    d(~mask) = 0;
end

d_flip = d*sin(tipangle/180*pi);
save d_flip d_flip
figure,im(d_flip);
axis off; colorbar; title('target excitation pattern');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set weighting matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [gbx gby gbz] = gradient(b0map);
% gb = gbx.^2+gby.^2 +gbz.^2;
gb = ones(size(b0map)); % just use all one here
weightIm = (gb).^0.5; % same dimension as image

if extype == 3
    % get the outer volume
    SEd = strel('disk',8);
    SEe = ones(3,3);
    %SEe = strel('disk',2);
    for i = 1:size(d,3)
        d_dilated(:,:,i) = imdilate(abs(d(:,:,i)),SEd);
        d_erode(:,:,i) = imerode(abs(d(:,:,i)),SEe);
    end
    iv_roi = (d_erode > 0.1);
    ov_roi = logical(1-(d_dilated > 0.1));
    edge_roi = logical(1 - iv_roi - ov_roi);
end
% weightIm(edge_roi) = 0.1*weightIm(edge_roi);
% weight = diag_sp(weightIm(roi));
% save weighting weight weightIm gb;
if 0
    ov_roi(:,:,z1-1) = ov_roi(:,:,z1);
    ov_roi(:,:,z2) = ov_roi(:,:,z2);
end

save rois iv_roi ov_roi edge_roi

weightIm(ov_roi) = 2.4; %2;  % very important to suppress outer-volume in steady-state imaging
weightIm(iv_roi) = 0.7; % 0.3;
weightIm(edge_roi) = 0.3; %0.02;

% weight = sparse(diag(weightIm(roi)));
weight = diag_sp(weightIm(roi));
save weighting weight weightIm gb;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  pulse design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [b1,gx,gy,gz] = spokerf(d,b0map,mask,tipangle,n_pe,slice_thickness,ishard,simuRange);
% [b1,gx,gy,gz] = ktpointsRF(d,b0map,mask,tipangle,n_pe,ishard,simuRange);
% [b1,gx,gy,gz] = ktContRF(d,b0map,mask,tipangle,n_pe,ishard,xyzRange);
% mask = ones(size(mask));
% [b1,gx,gy,gz,M_bloch] = spectralRF(d_flip,b0map,sens,roi,simuRange,12e-3);

%% compare with SPINS pulse
if trajtype == 1
    addpath ..
    addpath ../Gradient/
    spinfast = 1;
    if spinfast == 1
        % Nov2
        multi = 3;
        trajPara.dur   = 10e-3;       % Duration in seconds
        trajPara.kmax  = 50*multi;         % Maximum extent of k-space in rad m^-1
        trajPara.u     = 30*pi/trajPara.dur; %*multi;   % Polar angular velocity default 8*pi
        trajPara.v     = 20*pi/trajPara.dur; %*multi;   % Azimuthal angular velocity default 2*pi
        trajPara.alpha = 10;         % Speed of transition between slow and fast radial phase
        trajPara.beta  = 0.5;        % Position of transition between slow and fast radial phase

		% the following produces 4.0ms SPINS gradients
        trajPara.dur   = 9e-3;       % Duration in seconds
        trajPara.kmax  = 50*multi;         % Maximum extent of k-space in rad m^-1
        trajPara.u     = 25*pi/trajPara.dur; %*multi;   % Polar angular velocity default 8*pi
        trajPara.v     = 18*pi/trajPara.dur; %*multi;   % Azimuthal angular velocity default 2*pi
        trajPara.alpha = 10;         % Speed of transition between slow and fast radial phase
        trajPara.beta  = 0.5;        % Position of transition between slow and fast radial phase

			% JFN 
        trajPara.dur   = 9e-3/2;       % Duration in seconds
        trajPara.kmax  = 30*multi;         % Maximum extent of k-space in rad m^-1
        trajPara.u     = 25*pi/trajPara.dur; %*multi;   % Polar angular velocity default 8*pi
        trajPara.v     = 18*pi/trajPara.dur; %*multi;   % Azimuthal angular velocity default 2*pi
        trajPara.alpha = 10;         % Speed of transition between slow and fast radial phase
        trajPara.beta  = 0.5;        % Position of transition between slow and fast radial phase
    else
        multi = 1.5;
        trajPara.dur   = 4e-3;       % Duration in seconds
        trajPara.kmax  = 50*multi;         % Maximum extent of k-space in rad m^-1
        trajPara.u     = 8*pi/trajPara.dur; %*multi;   % Polar angular velocity default 8*pi
        trajPara.v     = 4*pi/trajPara.dur; %*multi;   % Azimuthal angular velocity default 2*pi
        trajPara.alpha = 10;         % Speed of transition between slow and fast radial phase
        trajPara.beta  = 0.5;        % Position of transition between slow and fast radial phase

        trajPara.dur   = 10e-3;       % Duration in seconds
        trajPara.kmax  = 50*multi;         % Maximum extent of k-space in rad m^-1
        trajPara.u     = 16*pi/trajPara.dur; %*multi;   % Polar angular velocity default 8*pi
        trajPara.v     = 8*pi/trajPara.dur; %*multi;   % Azimuthal angular velocity default 2*pi
        trajPara.alpha = 10;         % Speed of transition between slow and fast radial phase
        trajPara.beta  = 0.5;        % Position of transition between slow and fast radial phase
    end
    % d_flip = ones(size(d))*sin(tipangle/180*pi);
    d_flip = d*sin(tipangle/180*pi);
    %    [b1,gx,gy,gz,tmpIm,kspin] = spinsrf(d_flip,b0map,mask,sens,xyzRange,trajPara,0);    % spiral pulse
    [b1,gx,gy,gz] = spinsrf_fast(d_flip,b0map,mask,sens,xyzRange,trajPara,0);    % spiral pulse
	kspins = g2k_hao([gx(:) gy(:) gz(:)]);
	save kspins kspins
    %T1 = 1000; T2 = 100; %dt = 4e-6;
    M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt*1e-3,b0map_simu,mask,T1/1000,T2/1000);
    % figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default;
    M_spins = M_pr;
    nrmse_spin = norm(M_pr(roi)- d_flip(roi))/sqrt(sum(roi(:)))/max(abs(d_flip(roi)))
    diff_bloch = embed((M_pr(mask)-d_flip(mask)),mask);
    if extype ~=3
        figure,im(diff_bloch,[0 max(abs(d_flip(:)))],'cbar'),colormap default;
        title(sprintf('nrmse: %.3f', nrmse_spin));
    else
        figure,im(M_pr,[0 max(abs(d_flip(:)))],'cbar'),colormap default;
        title(sprintf('nrmse: %.3f', nrmse_spin));
    end
    tt = 1:length(b1);
    %dt = 4e-3;
    tt = tt*dt;
    length(b1)*dt
    
    % add crusher and write out to .wav
    ncycles = 8;
    zthick = 1;
    netarea = 0;
    gc = makecrusher_hao(ncycles, zthick, netarea); %1cm 8 cycles Do we have to include the netarea here?
    b1_nc = b1;
    gx_nc = gx;
    gy_nc = gy;
    gz_nc = gz;
    
    b1 = [zeros(length(gc),ncoils);b1];
    gx = [zeros(size(gc));gx];
    gy = [zeros(size(gc));gy];
    gz = [gc;gz];
    if spinfast == 1
        writewav_hao('spins_fast.wav',b1,[gx,gy,gz],tipangle);
        plotwav('spins_fast.wav');
    else
        writewav_hao('spins.wav',b1,[gx,gy,gz],tipangle);
        plotwav('spins.wav');
    end
end

%% compare with stack of spiral pulse
if trajtype == 2
    spiralPara.ncycle = 6; %dx = fov/ncycle;
    spiralPara.ncycleZ = 5;
    spiralPara.Fov = -xyzRange.x(1)*2;
    spiralPara.Fovz = 1*(xyzRange.z(end) - xyzRange.z(1) + xyzRange.z(2) - xyzRange.z(1));
    d_flip = d*sin(tipangle/180*pi);
    [b1,gx,gy,gz,projIm_SoS] = stack_spiral_rf(d_flip,b0map,mask,sens,xyzRange,spiralPara,weightIm);    % spiral pulse
    %T1 = 1000; T2 = 100; %dt = 4e-6;
    M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt*1e-3,b0map_simu,mask,T1/1000,T2/1000);
    % figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default;
    M_SoS = M_pr;
    nrmse_spin = norm(M_pr(roi)- d_flip(roi))/sqrt(sum(roi(:)))/max(abs(d_flip(roi)))
    diff_bloch = embed((M_pr(mask)-d_flip(mask)),mask);
    if extype ~=3
        figure,im(diff_bloch,[0 max(abs(d_flip(:)))],'cbar'),colormap default;
        title(sprintf('nrmse: %.3f', nrmse_spin));
    else
        figure,im(M_pr,[0 max(abs(d_flip(:)))],'cbar'),colormap default;
        title(sprintf('nrmse: %.3f', nrmse_spin));
    end
    tt = 1:length(b1);
    %dt = 4e-3;
    tt = tt*dt;
    length(b1)*dt
    
    % add crusher and write out to .wav
    ncycles = 8;
    zthick = 1;
    netarea = 0;
    gc = makecrusher_hao(ncycles, zthick, netarea); %1cm 8 cycles Do we have to include the netarea here?
    b1_nc = b1;
    gx_nc = gx;
    gy_nc = gy;
    gz_nc = gz;
	 ksos = g2k_hao([gx_nc(:) gy_nc(:) gz_nc(:)]);
    save ksos ksos
    
    b1 = [zeros(length(gc),ncoils);b1];
    gx = [zeros(size(gc));gx];
    gy = [zeros(size(gc));gy];
    gz = [gc;gz];
end

% design tipdown pulse
switch exptype
	case 'iv'
		[b1,gx,gy,gz,projIm,M_bloch,K_pro] = ktContRF_greedy2(d,b0map,b0map_simu,sens, mask,tipangle,n_pe,ishard,simuRange,tol,extype,isBssfp,reply,trajtype,niter);
	case 'slr'
		[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readmodv2('~jfnielse/projects/myelin/tipdown.mod');
		nomflip = paramsfloat(11);
		b1 = rho.*exp(1i*(theta-pi/2))/nomflip*tipangle;    % SLR pulse produces pi/2 phase, so remove that
end

save rfpulse b1 gx gy gz



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Done with pulse design.
%%  Display and write to file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(exptype,'iv')
	pulseLength = length(b1)*dt;
	max(abs(b1(:)))
	nrmse = norm(M_bloch(roi) - d_flip(roi))/sqrt(sum(roi(:)))/max(abs(d_flip(roi)))
	nrmse_proj = norm(projIm(roi)- d_flip(roi))/sqrt(sum(roi(:)))/max(abs(d_flip(roi)))
	% figure,im(projIm-d_flip,[0 max(abs(d_flip(:)))]); colormap default; title(sprintf('diff image; small-tip model; nrmse: %.3f', nrmse_proj));
	figure,im(projIm,[0 max(abs(d_flip(:)))]); colormap default; title(sprintf('small-tip model; nrmse: %.3f', nrmse_proj)); axis off;

	diff_bloch = embed((M_bloch(mask)-d_flip(mask)),mask);
	figure,im(diff_bloch,[0 max(abs(d_flip(:)))],'cbar'); colormap default;
	title(sprintf('error image; nrmse: %.3f', nrmse)); axis off;

	im_bloch = embed((M_bloch(mask)),mask);
	figure,im(im_bloch,[0 max(abs(d_flip(:)))],'cbar'); colormap default;
	title(sprintf('Bloch sim')); axis off;

	dtarget = embed(d_flip(mask),mask);
	dbloch = embed(M_bloch(mask),mask);
	save results dtarget dbloch

	if extype ~=3
  	  figure,im('notick',diff_bloch,[0 max(abs(d_flip(:)))],'cbar',sprintf('error image; nrmse: %.3f', nrmse)),colormap default;
	else
		figure,im('notick',M_bloch,[0 max(abs(d_flip(:)))],'cbar',sprintf('nrmse: %.3f', nrmse)),colormap default;
		max_ov = max(abs(M_bloch(ov_roi)))/max(abs(d_flip(roi)))
	end
end

if strcmp(exptype,'iv')
    if simumode == 2
		fname = sprintf('tipdown.mod');
        if isBssfp
            mat2wav(abs(b1),angle(b1),gx,gy,gz,tipangle,fname,'balanced tipdown');
        else
            mat2wav(abs(b1),angle(b1),gx,gy,gz,tipangle,'propose_sos3d.wav','tipdown');
            plotwav_sos3d('propose_sos3d.wav');
        end
    end
end
