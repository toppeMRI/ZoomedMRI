
% main file for kt-points
close all;
clear;
setupSimu_ktpoints; 

slices = 70; 
subjects = [1];
n_slices = length(subjects)*length(slices);
ex_type = 1; % 1: uniform; 2: prephaseing; 3: rectangle
miniMax = 0; % In Daehyun's method, the structure of algo. is different. So the minimax here is not what I proposed, but replace all the L2 minimization in Daehyun's code with Linf minimization 
n_pe = 20;
ncand = 10; 
ishard = 0; 
type = 2; % calculate correlation: 1:omp-l1, 2:omp-l2
lambda = 20; %regulation term


masks = zeros(64,64,n_slices);
fieldmaps = zeros(64,64,n_slices);
densities = zeros(64,64,n_slices);
indices = 1;
i_sl = 1; 
i_sb = 1; 
[d,sens,mask,density,dx,dy,sx,sy,ncoils,fieldmap,nfreq,pulse_params] = loadup_simulation(ex_type, slices(i_sl), subjects(i_sb));
d = zeros(size(d)); 
d(25:40,25:40) = 1; 
masks(:,:,indices)=mask;
fieldmaps(:,:,indices)=fieldmap;
densities(:,:,indices)=density;

figure,im(masks),colormap gray; title('mask');
figure,im(densities),colormap gray; title('density');
figure,im(fieldmaps),colormap default; title('fieldmap');




% our proposed method
kx_proposed = zeros(n_pe,n_slices); 
ky_proposed = zeros(n_pe,n_slices);
times_proposed = zeros(n_slices,1);

sl_index = 1;


% normalized sensitivity to determine phase encoding locations
sens_norm = zeros(size(sens));

% % normalize sensitivity for determining phase encoding locations
% for i = 1:ncoils
%     sens_norm(:,:,i)=sens(:,:,i);
%     sens_i = squeeze(sens(:,:,i));
%     sens_norm(:,:,i) = sens_norm(:,:,i)/norm(sens_i(:));
% end
sens_norm = ones(64,64); !!!
sens = sens_norm;  !!!

[nrmse,kx_proposed,ky_proposed] = determine_pe_field_greedy3(d,sens_norm,mask,n_pe,type,ncand,fieldmap,sx,sy,pulse_params, miniMax,ishard);


%% plot proposed method: 
sl_index = 1;

x_range = ([0:1:dx-1]-floor(dx/2))*sx/dx;
y_range = ([0:1:dy-1]-floor(dy/2))*sy/dy;
z_range = [0];


kx = kx_proposed(:,sl_index);
ky = ky_proposed(:,sl_index);
[x_est,projection,pr_nrmse] = determine_weights_field2(sx,sy,kx,ky,d,mask,sens,fieldmap,pulse_params,lambda, miniMax);        


M_pr = zeros(dy,dx,length(z_range));

[b1,gx,gy,gz] = create_pulse(x_est,n_pe,ncoils,kx,ky,pulse_params,dx,dy,sx,sy,ishard); %output is in Tesla
disp(sprintf('the net pulse length %f msec',size(b1,1)*pulse_params.dt*1000));

save pulsesoutput; 
for i = 1:length(z_range)
    M_pr(:,:,i)=parallel_blochsim_field(0,b1,gx,gy,gz,sens,x_range,y_range,z_range(i),pulse_params.dt,fieldmap,0);
end

M_pr_tr = sum(M_pr,3);

figure,im(abs(d)); colormap default; title('target phase pattern'); 
figure,im(abs(projection)*sin(pulse_params.tip_angle*2*pi/360));colormap default; title('projection'); 
figure,im(abs(M_pr_tr));colormap default; title('excitation pattern (Bloch)'); 

tip_angle = pulse_params.tip_angle; 
d_tar = sin(tip_angle*pi/180); 
diffIm = M_pr_tr - d_tar;
diffIm(~mask) = 0; 
%figure,im(diffIm); 
max(abs(diffIm(:)))./d_tar; 

tt = 1:length(b1);
dt = 4e-3; 
tt = tt*dt; 
% figure,subplot(4,1,1),plot(tt, abs(b1(:,1))*1e4); ylabel('gauss'); 
% subplot(4,1,2),plot(tt,abs(gx)); ylabel('g/cm'); 
% subplot(4,1,3),plot(tt,abs(gy)); 
% subplot(4,1,4),plot(tt,abs(gz)); ylim([-5,5]); 
% 




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
