close all;
clear;
%addpath ../
%addpath ../split
%slices = [25:10:45];
%subjects = [1:2];
% slices = [25:10:35];
slices = 50; 
subjects = [1];
n_slices = length(subjects)*length(slices);
ex_type = 2; % 1: uniform; 2: prephaseing; 3: rectangle
miniMax = 0; % In Daehyun's method, the structure of algo. is different. So the minimax here is not what I proposed, but replace all the L2 minimization in Daehyun's code with Linf minimization 
n_pe = 10;
ncand = 10; 
%% Figure 1

load sensitivities;
% sens = cat(3, sens(:,:,1), sens(:,:,2)); % Hao
figure,im(abs(sens)),colormap gray;
title('figure 1, sensitivity magnitude');


%% Figure 2
%%{T
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
%fieldmap = 0*fieldmap; !!!!
figure,im(masks),colormap gray; title('mask'); 
figure,im(densities),colormap gray; title('density'); 
figure,im(fieldmaps),colormap default; title('fieldmap'); 


%% PE selection for all methods
%%{
multiple_rate = 0.01;
error_bound = 10^(-9); 
type = 2; %1:omp-l1, 2:omp-l2
lambda = 0; %regulation term

% our proposed method
kx_proposed = zeros(n_pe,n_slices); 
ky_proposed = zeros(n_pe,n_slices);
times_proposed = zeros(n_slices,1);

sl_index = 1;
for i_sb = 1:length(subjects)
    for i_sl = 1:length(slices)

        [d,sens,mask,density,dx,dy,sx,sy,ncoils,fieldmap,nfreq,pulse_params] = loadup_simulation(ex_type, slices(i_sl), subjects(i_sb));

        % normalized sensitivity to determine phase encoding locations
        sens_norm = zeros(size(sens));

        % normalize sensitivity for determining phase encoding locations
        for i = 1:ncoils
            sens_norm(:,:,i)=sens(:,:,i);
            sens_i = squeeze(sens(:,:,i));
            sens_norm(:,:,i) = sens_norm(:,:,i)/norm(sens_i(:));
        end
        sens_norm = ones(64,64); !!!
        sens = sens_norm;  !!!
        disp(sprintf('slice %d out of %d',sl_index,n_slices));        
        disp('proposed method');
        tic;
        [nrmse,kx_proposed(:,sl_index),ky_proposed(:,sl_index)] = determine_pe_field_greedy3(d,sens_norm,mask,n_pe,type,ncand,fieldmap,sx,sy,pulse_params, miniMax);
        times_proposed(sl_index) = toc;

        kx = kx_proposed(:,sl_index);
        ky = ky_proposed(:,sl_index);
        [x_est,projection,pp_nrmse] = determine_weights_field2(sx,sy,kx,ky,d,mask,sens,fieldmap,pulse_params,lambda, miniMax);        
        pp_nrmse
        
        sl_index = sl_index + 1;
    end
end

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

save pulsesoutput; 
for i = 1:length(z_range)
    M_pr(:,:,i)=parallel_blochsim_field(0,b1,gx,gy,gz,sens,x_range,y_range,z_range(i),pulse_params.dt,fieldmap,0);
end

M_pr_tr = sum(M_pr,3);

figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default; title('projection'); 
figure,im(abs(M_pr_tr));colormap default; title('excitation pattern (Bloch)'); 

tip_angle = pulse_params.tip_angle; 
d_tar = sin(tip_angle*pi/180); 
diffIm = M_pr_tr - d_tar;
diffIm(~mask) = 0; 
%figure,im(diffIm); 
max(abs(diffIm(:)))./d_tar

%figure,plot(abs(b1(:,1)));
%figure,plot(abs(gx)); 
%figure,plot(abs(gy)); 
%figure,plot(abs(gz)); ylim([-5,5]); 





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
