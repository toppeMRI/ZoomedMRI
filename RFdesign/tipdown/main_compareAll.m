% main file for continuous trajectory based on kt points
%% setup
clear;
rFOV = 1;

for doOpt = 0:1
    for iniTraj = 1:4
        close all;
        setupSimu_Cont; !!! check the b0map is scaled or not
        % b0map = 1e-10*b0map; %  mainly for testing
        % b0map_simu = 1e-10*b0map_simu;
        % save b0map b0map;
        % save b0map_simu b0map_simu;
        
        dospin = 0; % compare with spins trajectory
        doSoS = 0; % compare with stack of spiral trajectory
        tipangle = 10;
        isBssfp = 0; % balanced gradient
        if isBssfp
            tipangle = 2*tipangle;
        end
        ishard = 1; % used to determine the subpulse type (spoke only) and its duration (only for the b0 correction)
        tol = 0.5;
        
      
        extype = 3; % 1: uniform excitation; 2: pre-phaseing; 3: cube/hocky pot; 4: slice selective pre-phasing;
        if rFOV == 0
            extype = 2; 
        end 
        
        dorect = 1;
        dosmooth = 0; % smooth the edge of target excitation pattern
        
        if extype == 1
            d = ones(size(b0map_simu));
            d(~mask) = 0;
        end
        
        if extype == 2 % pre-phaseing;
            n_pe = 60;
            %    n_pe = 27; % for kt-points (to achieve the same pulse length);
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
        
        
        if extype == 3
            % n_pe = 65; % nov 2 phantom
            n_pe = 70; % for ktraj paper
            %     n_pe = 21; % for 4.5ms ktpoints
            %    n_pe = 100;
            d = zeros(size(b0map_simu));
            if dorect
                %       d(25:40,25:40,ceil(size(b0map_simu,3)/2):ceil(size(b0map_simu,3)/2)+2) = 1; % nov2
                d(25:40,25:40,ceil(size(b0map_simu,3)/2):ceil(size(b0map_simu,3)/2)+2) = 1;
                %          d(20:35,20:45,ceil(size(b0map_simu,3)/2-1):ceil(size(b0map_simu,3)/2)+2) = 1; % test
            else
                xr = 1:size(d,1);
                yr = 1:size(d,2);
                [xxtmp,yytmp] = ndgrid(xr-32,yr-32);
                radius = 4.5/24*size(d,1);
                circInd = xxtmp.^2 + yytmp.^2 < radius.^2;
                nz = size(b0map_simu,3);
                z1 = ceil(nz/2)-1;
                z2 = ceil(nz/2)+1;
                circInd = repmat(circInd, [1 1 z2-z1+1]);
                d(:,:,z1:z2) = circInd;
                if dosmooth == 1
                    for i = 1:size(d,3)
                        h = fspecial('gaussian', 5, 1);
                        d(:,:,i) = imfilter(d(:,:,i), h);
                    end
                end
            end
        end
        
        if extype == 4 %slice selective pre-phasing;
            n_pe = 60;
            dp = exp(1i*b0map_simu*Tread/1000*2*pi);
            d = zeros(size(dp));
            d(:,:,ceil(size(b0map_simu,3)/2-1):ceil(size(b0map_simu,3)/2)+3) = dp(:,:,ceil(size(b0map_simu,3)/2-1):ceil(size(b0map_simu,3)/2)+3);
            d(~mask) = 0;
        end
        
        
        d_flip = d*sin(tipangle/180*pi);
        figure,im(d_flip);
        axis off; colorbar; title('target excitation pattern');
        
        % setting weighting matrix
        % [gbx gby gbz] = gradient(b0map);
        % gb = gbx.^2+gby.^2 +gbz.^2;
        gb = ones(size(b0map)); % just use all one here
        weightIm = (gb).^0.5; % same dimension as image
        if extype == 3
            weightIm = 1*(1 - d)+d;
            % figure,im(weightIm);
        end
        % weight = sparse(diag(weightIm(roi)));
        weight = diag_sp(weightIm(roi));
        save weighting weight weightIm gb;
        
        if extype == 3
            % get the outer volume
            SEd = strel('disk',8);
            SEe = ones(3,3);
            for i = 1:size(d,3)
                d_dilated(:,:,i) = imdilate(d(:,:,i),SEd);
                d_erode(:,:,i) = imerode(d(:,:,i),SEe);
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
        
        %%  pulse design
        % [b1,gx,gy,gz] = spokerf(d,b0map,mask,tipangle,n_pe,slice_thickness,ishard,simuRange);
        % [b1,gx,gy,gz] = ktpointsRF(d,b0map,mask,tipangle,n_pe,ishard,simuRange);
        % [b1,gx,gy,gz] = ktContRF(d,b0map,mask,tipangle,n_pe,ishard,xyzRange);
        % mask = ones(size(mask));
        
        
        % [b1,gx,gy,gz,M_bloch] = spectralRF(d_flip,b0map,sens,roi,simuRange,12e-3);
        
        [b1,gx,gy,gz,projIm,M_bloch,K_pro] = kb_compareAll(d,b0map,sens, mask,tipangle,n_pe,ishard,simuRange,tol,extype,isBssfp,doOpt,iniTraj,rFOV);
        length(b1)*dt
        max(abs(b1(:)))
        nrmse = norm(M_bloch(roi) - d_flip(roi))/sqrt(sum(roi(:)))/max(abs(d_flip(roi)))
        nrmse_proj = norm(projIm(roi)- d_flip(roi))/sqrt(sum(roi(:)))/max(abs(d_flip(roi)))
        % figure,im(projIm-d_flip,[0 max(abs(d_flip(:)))]); colormap default; title(sprintf('diff image; small-tip model; nrmse: %.3f', nrmse_proj));
        figure,im(projIm,[0 max(abs(d_flip(:)))]); colormap default; title(sprintf('small-tip model; nrmse: %.3f', nrmse_proj)); axis off;
        
        diff_bloch = embed((M_bloch(mask)-d_flip(mask)),mask);
        figure,im(diff_bloch,[0 max(abs(d_flip(:)))],'cbar'); colormap default;
        title(sprintf('error image; nrmse: %.3f', nrmse)); axis off;
        
        if extype ~=3
            figure,im('notick',diff_bloch,[0 max(abs(d_flip(:)))],'cbar',sprintf('error image; nrmse: %.3f', nrmse)),colormap default;
            %    ylabel('SoS');
            %    xlabel('Initialization');
            %    xlabel('After interior-point optimization');
            %title(sprintf('bloch simulation nrmse: %.3f', nrmse)); axis off;
        else
            figure,im('notick',M_bloch,[0 max(abs(d_flip(:)))],'cbar',sprintf('nrmse: %.3f', nrmse)),colormap default;
            %    title(sprintf('nrmse: %.3f', nrmse)); axis off;
            max_ov = max(abs(M_bloch(ov_roi)))/max(abs(d_flip(roi)))
        end
    
    if rFOV == 1
        save(['forPaper_rFOV_' 'doOpt' num2str(doOpt) '_ktraj' num2str(iniTraj)], 'projIm', 'd_flip', 'roi', 'b0map', 'M_bloch', 'diff_bloch', 'nrmse','nrmse_proj','b1', 'gx', 'gy', 'gz'); 
    else
        save(['forPaper_prephasing_' 'doOpt' num2str(doOpt) '_ktraj' num2str(iniTraj)], 'projIm', 'd_flip', 'roi', 'b0map', 'M_bloch', 'diff_bloch', 'nrmse','nrmse_proj','b1', 'gx', 'gy', 'gz'); 
    end
    end
end

%% add crusher and write out to .wav
if isBssfp == 0
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
end
% if 0 % write wavs for stfr5.e
%     if simumode == 2
%         writewav_hao('propose.wav',b1,[gx,gy,gz],tipangle);
%         plotwav('propose.wav');
%     end
% end
%
% if 0 % generate gradient wavs for stfr5.e
%     writewav_hao('propose_gx.wav',b1,[gx,zeros(size(gy)),zeros(size(gz))],tipangle);
%     writewav_hao('propose_gy.wav',b1,[gy,zeros(size(gy)),zeros(size(gz))],tipangle);
%     gz(1:length(gc)) = 0;
%     writewav_hao('propose_gz.wav',b1,[gz,zeros(size(gy)),zeros(size(gz))],tipangle);
%     writewav_hao('propose_g_ref.wav',b1,[zeros(size(gx)),zeros(size(gy)),zeros(size(gz))],tipangle);
% end

if 1 % write wavs for stfr sos3d.e
    if simumode == 2
        if isBssfp
            mat2wav(abs(b1),angle(b1),gx,gy,gz,tipangle,'propose_bssfp.wav','balanced tipdown');
            plotwav_sos3d('propose_bssfp.wav');
            
        else
            mat2wav(abs(b1),angle(b1),gx,gy,gz,tipangle,'propose_sos3d.wav','tipdown');
            %         writewav_sos3d('propose_sos3d.wav',b1,[gx,gy,gz],tipangle);
            plotwav_sos3d('propose_sos3d.wav');
            
        end
    end
end

if 0 % generate gradient wavs for stfr5.e
    writewav_sos3d('propose_gx.wav',b1,[gx,zeros(size(gy)),zeros(size(gz))],tipangle);
    writewav_sos3d('propose_gy.wav',b1,[gy,zeros(size(gy)),zeros(size(gz))],tipangle);
    gz(1:length(gc)) = 0;
    writewav_sos3d('propose_gz.wav',b1,[gz,zeros(size(gy)),zeros(size(gz))],tipangle);
    writewav_sos3d('propose_g_ref.wav',b1,[zeros(size(gx)),zeros(size(gy)),zeros(size(gz))],tipangle);
end

%% load wav in Nov2,2013 and do simulation
if 0
    clear;
    setupSimu_Cont;
    d = zeros(size(b0map_simu));
    d(25:40,25:40,ceil(size(b0map_simu,3)/2):ceil(size(b0map_simu,3)/2)+2) = 1;
    %    d(25:40,25:40,ceil(size(b0map_simu,3)/2):ceil(size(b0map_simu,3)/2)+4) = 1;
    
    d_flip = d*sin(10/180*pi);
    xr = 10:55;
    yr = 10:55;
    zr = 9:2:20;
    d_flip = d_flip(:,:,5:end-2);
    resoz = FOVz/(Nz); %
    xyzRangeAll.z = (-Nz/2:Nz/2-1)*resoz; % here simulation is set to match with pulse design (cont_gen_rf.m) before I fix the 0.5 shift problem.
    scale = 0.9;
    simuRange.x = simuRange.x/scale;
    simuRange.y = simuRange.y/scale;
    %    resoz = FOVz/(Nz-1); %
    %    xyzRangeAll.z = -FOVz/2:resoz:FOVz/2;
    simuRange.z = xyzRangeAll.z(zr);
    sens = ones(1,64,64,length(zr));
    b0map = fieldmap(:,:,zr);
    mask = roiAll(:,:,zr);
    
    %  a interesting thing: SPINS traj do not need to negate gx, but the proposed
    %  need.
    %       [desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav('spins_fast.wav');
    %       b1 = rho.*exp(1i*theta);
    %       M_pr = parallel_blochCim(0,conj(b1),gx,-gy,-gz,sens,simuRange.x,simuRange.y,simuRange.z,4e-6,b0map,mask,1,0.1);
    [desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav('proposeNov2_2.wav');
    %        [desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav('proposeNov2GD2.wav');
    %     [desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav('ktpoints.wav'); % May 03, 14
    b1 = rho.*exp(1i*theta);
    M_pr = parallel_blochCim(0,conj(b1),-gx,gy,gz,sens,simuRange.x,simuRange.y,simuRange.z,4e-6,b0map,mask,1,0.1);
    % scanner has negative gx.
    nrmse = norm(M_pr(:)- d_flip(:))/(sqrt(sum(double(mask(:)))))/ sin(10/180*pi)
    
    figure,im('row',2,M_pr(xr,yr,:),[0, sin(10/180*pi)], 'cbar'), colormap default; axis off; title('');
end
%% evaluate the effect of k traj deviation
if 0 % compare the simulation result using real k traj.
    load measuredKG;
    gamma = 4257;
    dt = 4e-6;
    g = [gx_nom,gy_nom,gz_nom];
    k = -flipud(cumsum(flipud(g))*dt*gamma);
    [b1,gx1,gy1,gz1,projIm] = cont_gen_rf(d_flip,b0map,mask,sens,simuRange,k,g);
    T1 = 1000; T2 = 100; dt = 4e-6;
    M_pr = parallel_blochCim(0,b1,gx1,gy1,gz1,sens,simuRange.x,simuRange.y,simuRange.z,dt,b0map,mask,T1/1000,T2/1000);
    nrmse = norm(M_pr(roi)- d_flip(roi))/sqrt(sum(roi(:)))/max(d_flip(roi))
    figure,im(M_pr); title(sprintf('nrmse: % .3g',nrmse));
end


%% compare with SPINS pulse
if dospin == 1
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
    else
        multi = 1.5;
        trajPara.dur   = 4e-3;       % Duration in seconds
        trajPara.kmax  = 50*multi;         % Maximum extent of k-space in rad m^-1
        trajPara.u     = 8*pi/trajPara.dur; %*multi;   % Polar angular velocity default 8*pi
        trajPara.v     = 4*pi/trajPara.dur; %*multi;   % Azimuthal angular velocity default 2*pi
        trajPara.alpha = 10;         % Speed of transition between slow and fast radial phase
        trajPara.beta  = 0.5;        % Position of transition between slow and fast radial phase
    end
    % d_flip = ones(size(d))*sin(tipangle/180*pi);
    d_flip = d*sin(tipangle/180*pi);
    %    [b1,gx,gy,gz,tmpIm,kspin] = spinsrf(d_flip,b0map,mask,sens,xyzRange,trajPara,0);    % spiral pulse
    [b1,gx,gy,gz] = spinsrf_fast(d_flip,b0map,mask,sens,xyzRange,trajPara,0);    % spiral pulse
    T1 = 1000; T2 = 100; dt = 4e-6;
    M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt,b0map_simu,mask,T1/1000,T2/1000);
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
    dt = 4e-3;
    tt = tt*dt;
    % figure,subplot(4,1,1),plot(tt, abs(b1(:,1))); ylabel('gauss');
    % subplot(4,1,2),plot(tt,abs(gx)); ylabel('g/cm');
    % subplot(4,1,3),plot(tt,abs(gy));
    % subplot(4,1,4),plot(tt,abs(gz)); xlabel('ms');
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
if doSoS == 1
    spiralPara.ncycle = 6; %dx = fov/ncycle;
    spiralPara.ncycleZ = 5;
    spiralPara.Fov = -xyzRange.x(1)*2;
    spiralPara.Fovz = 1*(xyzRange.z(end) - xyzRange.z(1) + xyzRange.z(2) - xyzRange.z(1));
    d_flip = d*sin(tipangle/180*pi);
    [b1,gx,gy,gz,projIm_SoS] = stack_spiral_rf(d_flip,b0map,mask,sens,xyzRange,spiralPara,weightIm);    % spiral pulse
    T1 = 1000; T2 = 100; dt = 4e-6;
    M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt,b0map_simu,mask,T1/1000,T2/1000);
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
    dt = 4e-3;
    tt = tt*dt;
    % figure,subplot(4,1,1),plot(tt, abs(b1(:,1))); ylabel('gauss');
    % subplot(4,1,2),plot(tt,abs(gx)); ylabel('g/cm');
    % subplot(4,1,3),plot(tt,abs(gy));
    % subplot(4,1,4),plot(tt,abs(gz)); xlabel('ms');
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

%% plot figures for the proposal
if 0
    figure,im(M_bloch(:,:,7:end-6),[0 max(abs(d_flip(:)))],'cbar'),colormap default;
    title(sprintf('nrmse: %.3f', nrmse));
    print -depsc proposal_propose_2ms
    figure,im(M_spins(:,:,7:end-6),[0 max(abs(d_flip(:)))],'cbar'),colormap default;
    title(sprintf('nrmse: %.3f', nrmse_spin));
    print -depsc proposal_spins_26ms
    figure,im(M_pr(:,:,7:end-6),[0 max(abs(d_flip(:)))],'cbar'),colormap default;
    title(sprintf('nrmse: %.3f', nrmse_sos));
    print -depsc proposal_sos_28ms
end