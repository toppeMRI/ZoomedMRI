close all;
clear;
%slices = [25:10:45];
%subjects = [1:2];
% slices = [25:10:35];
slices = 50;
subjects = [1];
n_slices = length(subjects)*length(slices);
ex_type = 1;
miniMax = 0;

%% Figure 1

% load sensitivities;
% % sens = cat(3, sens(:,:,1), sens(:,:,2)); % Hao
% figure,im(abs(sens)),colormap gray;
% title('figure 1, sensitivity magnitude');


%% Figure 2
%{T
masks = zeros(64,64,n_slices);
fieldmaps = zeros(64,64,n_slices);
densities = zeros(64,64,n_slices);
indices = 1;
for i_sb = 1:length(subjects)
    for i_sl = 1:length(slices)
        [d,sens,mask,density,dx,dy,sx,sy,ncoils,fieldmap,nfreq,pulse_params] = loadup_simulation(ex_type, slices(i_sl), subjects(i_sb));
        
        masks(:,:,indices)=mask;
        fieldmaps(:,:,indices)=fieldmap;
        densities(:,:,indices)=density;
        indices = indices + 1;
    end
end

figure,im(masks),colormap gray;
figure,im(densities),colormap gray;
figure,im(fieldmaps),colormap default;


%% simulate all methods
%%{
n_peArray = [1:10];
for method = 1:3
    if method == 1
        miniMax = 0;
        ncand = 1; 
    elseif method == 2
        miniMax = 0;
        ncand = 1; 
    else 
        miniMax = 1; 
        ncand = 1; 
    end
    for ipe = 1:length(n_peArray)
        n_pe = n_peArray(ipe);
        multiple_rate = 0.01;
        error_bound = 10^(-9);
        type = 2; %1:omp-l1, 2:omp-l2
        lambda = 0; %regulation term
        % n_pe = 3;
        % our proposed method
        kx_proposed = zeros(n_pe,n_slices);
        ky_proposed = zeros(n_pe,n_slices);
        times_proposed = zeros(n_slices,1);
        
        sl_index = 1;
        i_sl = 1;
        subjects = 1;
        [d,sens,mask,density,dx,dy,sx,sy,ncoils,fieldmap,nfreq,pulse_params] = loadup_simulation(ex_type, slices(i_sl), subjects(i_sb));
        
        disp(sprintf('slice %d out of %d',sl_index,n_slices));
        disp('proposed method');
        tic;
        [nrmse,kx_proposed(:,sl_index),ky_proposed(:,sl_index)] = determine_pe_field_greedy3_hao(d,sens,mask,n_pe,type,ncand,fieldmap,sx,sy,pulse_params, miniMax);
        times_proposed(sl_index) = toc;
        
        kx = kx_proposed(:,sl_index);
        ky = ky_proposed(:,sl_index);
        [x_est,projection,pp_nrmse] = determine_weights_field2(sx,sy,kx,ky,d,mask,sens,fieldmap,pulse_params,lambda, miniMax);
        
        
        %% plot proposed method:
        sl_index = 1;
        lambda = 0;
        [d,sens,mask,density,dx,dy,sx,sy,ncoils,fieldmap,nfreq,pulse_params] = loadup_simulation(ex_type, slices(1), 2);
        
        x_range = ([0:1:dx-1]-floor(dx/2))*sx/dx;
        y_range = ([0:1:dy-1]-floor(dy/2))*sy/dy;
        z_range = [0];
        
        
        kx = kx_proposed(:,sl_index);
        ky = ky_proposed(:,sl_index);
        [x_est,projection,pr_nrmse] = determine_weights_field2(sx,sy,kx,ky,d,mask,sens,fieldmap,pulse_params,lambda, miniMax);
        
        
        M_pr = zeros(dy,dx,length(z_range));
        
        [b1,gx,gy,gz] = create_pulse(x_est,n_pe,ncoils,kx,ky,pulse_params,dx,dy,sx,sy);
        disp(sprintf('the net pulse length %f msec',size(b1,1)*pulse_params.dt*1000));
        
        for i = 1:length(z_range)
            M_pr(:,:,i)=parallel_blochsim_field(0,b1,gx,gy,gz,sens,x_range,y_range,z_range(i),pulse_params.dt,fieldmap,0);
        end
        
        M_pr_tr = sum(M_pr,3);
        
        tip_angle = pulse_params.tip_angle;
        d_tar = sin(tip_angle*pi/180);
        diffIm = projection*d_tar - d_tar;
        %diffIm = M_pr_tr - d_tar;
        diffIm(~mask) = 0;
        
        maxErr_all(ipe,method) = max(abs(diffIm(:)))./d_tar;
        
        
    end
end
figure(31),plot(maxErr_all(:,1),'g'); hold on; plot(maxErr_all(:,2),'b'); hold on; plot(maxErr_all(:,3),'r'); 
% figure,im(diffIm);
% figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default;
% figure,im(abs(M_pr_tr));colormap default;
% figure,im(abs(M_pr_tr));colormap gray;




% kx_range = [-5 5]/sx;
% ky_range = [-5 5]/sy;
% figure;
% hold on;
% plot(kx/sx,ky/sy,'o-','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',11)
% plot(kx(1)/sx,ky(1)/sy,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',11);
% plot(kx(end)/sx,ky(end)/sy,'o','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',11);
% hold off;
% grid on;
% xlabel('kx (cycle/cm)');
% ylabel('ky (cycle/cm)');
% xlim(kx_range);
% ylim(ky_range);
% title(sprintf('Proposed method : selected PE locations and their ordering'));



return;
