function [B1,gx,gy,gz,projIm] = cont_gen_rf(d,fmap,roi,sens,xyzRange,contK, contG, B1_init, W)
% Based on: (1) Spiral RF pulse design, small-tip (Grissom MRM2005)
%           (2) Generate RF waveforms given k-traj
%
% Units:
% k-traj: cycle/cm
% gradient: Gauss/cm
% B1: Gauss
%
% Hao: Jun.25. 2012 Part of the code is borrow from spiralrf.m (Grissom)


%tipangle = 16; % pulse flip angle

Nz = size(fmap,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('sens')
   %load fdtdsens
   if 0
      load sens;  % from ~/projects/parallelrf/data/Apr2011/9/8coil_sens.mat
      %mask = sum(abs(sens),3) > 0.1*max(col(sum(abs(sens),3))) & roi;
   else
      sens = ones(1,64,64,Nz);
   end
end

% Hao: reshape the coil dimension as the 4th dimension
sens_tmp = [];
if size(sens,1) < 9
   for ic = 1:size(sens,1)
      sens_tmp(:,:,:,ic) = sens(ic,:,:,:);
   end
elseif length(size(sens)) < 4
   sens_tmp = ones(64,64,Nz,1);
end
sens = sens_tmp;

%% load mask;
%mask = logical(ones(size(roi)));
%sens(:,:,2:2:end) = [];
%save mask mask
load mask;
load roi;
dim = length(xyzRange.x); % dimension of square x-y grid
dimz = length(xyzRange.z);
Nc = size(sens,4); % number of Tx coils
FOV = -xyzRange.x(1)*2;
FOVz = -xyzRange.z(1)*2;
FOVz = xyzRange.z(end) - xyzRange.z(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the desired flip angle pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xx, yy, zz] = ndgrid(xyzRange.x,xyzRange.y,xyzRange.z);

sens = sens/max(abs(sens(:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a k-space trajectory and its gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = contK;
g = contG;


%g = g/10; %change from mt/m to Gauss/cm
Nt = size(g,1); % number of samples in pulses
NN = [Nt - size(k,1), size(k,1)];

dt = 4e-6; % seconds, RF and gradient sampling period



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the small-tip object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('Gsml')
   % field-map params
   tt = 0:dt:(Nt-1)*dt;tt = tt-(Nt-1)*dt;L = 4; %fmap = []; % Hz
   % nufft params
   J = 6;K = 2*dim;Kz = 2*dimz; Jz = 6;
   %nufft_args = {[dim dim],[J J],[K K],[dim dim]/2,'minmax:kb'};
   nufft_args = {[dim dim dimz],[J J Jz],[K K Kz],[dim dim dimz-1]/2,'minmax:kb'};
   
   gambar = 4257.6;	% gamma/2pi in Hz/g
   gam = gambar*2*pi; % gamma in radians/g
   % trick: the system matrix is just the transpose of a SENSE image recon matrix!
   %Gsml = Gmri_SENSE(k,mask,'fov',[FOV FOV],'basis',{'dirac'}, ...
   
   save tmpCont
   load tmpCont

   sens_reshape = reshape(sens,[dim*dim*dimz Nc]);
   sens_masked = sens_reshape(roi,:);
   FOVz_s = FOVz*dimz/(dimz-1); % nufft doesn't take fov and dim instead of xyzRange as input, so scale and shift to match with simulation
   %    FOVz_s = FOVz;
   Gsml = Gmri_SENSE(k,logical(mask),'fov',[FOV FOV FOVz_s],'basis',{'dirac'}, ...
      'nufft',nufft_args,'exact',0,... % try exact = 1 here;
      'sens',conj(sens_masked)*(-1i*gam*dt), ...
      'ti',tt,'L',L,'zmap',2*pi*1i*fmap)';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design the initial pulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('B1')
   % get the tikhonov penalty vector
   beta = 8.0; %change
   %     beta = 0;
   betavec = ones(Nc*Nt,1)*sqrt(beta);
   % penalize RF during the rewinder to force it to zero
   if 0  % forward spiral
      betavec = betavec+1000000*kron(ones(Nc,1),[zeros(NN(1),1);ones(NN(2),1)]);
   else
      [size(betavec) size(kron(ones(Nc,1),[zeros(NN(2),1)]))]
      betavec = betavec+1000000*kron(ones(Nc,1),[zeros(NN(2),1)]);
   end
   
   betavec = betavec+1000000*kron(ones(Nc,1),[ones(NN(1),1);zeros(NN(2),1)]);
   
   R = diag_sp(betavec);
   
   niters = floor(Nt/8);
   %niters = 20;
   disp 'Designing initial pulse'
   %[xS,info] = qpwls_pcg(zeros(Nc*Nt,1),Gsml,1,d*tipangle*pi/180,0, ...
   %		R,1,niters); %,mask);
   %		R,1,niters,ones(size(d)));
   if 1
      if ~isvar('B1_init') || B1_init == 0
         B1_init = zeros(Nc*Nt,1);
      end
      %[xS,info] = qpwls_pcg1(zeros(Nc*Nt,1),Gsml,1,d(mask)*tipangle*pi/180, ...
      
      %         if ~isvar(W)
      if 0
         W = diag_sp(ones(sum(roi(:)),1));
      else
         load weighting;
         W = weight;
         %             figure, im(weightIm); colormap default; title('weighting');
      end
      disp('running least squares');
      
      
      [xS,info] = qpwls_pcg1_hao(B1_init(:),Gsml,W,d(mask), ...
         R,'niter',niters);
   else
      [xS] = cgnr_jfn(Gsml,d(mask)*tipangle*pi/180,zeros(Nc*Nt,1),niters);   % for testing
   end
   
   B1(:,1) = xS(:,end);
end
size(B1)

% figure,plot(info.cost);
projIm = embed(Gsml*B1,mask); % see the image based on small tip model
% figure,im(projIm,'cbar'), colormap default;


% simulate


B1_init = reshape(B1(:,1),[Nt Nc]);
B1 = reshape(B1(:,1),[Nt Nc]);

gx = g(:,1);
gy = g(:,2);
gz = g(:,3);
B1 = [B1; zeros(length(gx)-length(B1),Nc);];

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% plot channel 1's pulses
subplot(221)
plot(0:dt*1000:(NN(1)-1)*dt*1000,abs(B1_init(1:NN(1),1)));
hold on
%plot(0:dt*1000:(NN(1)-1)*dt*1000,abs(B1_final(1:NN(1),1)),'r');
xlabel 'Time (ms)'
ylabel '|b_1(t)| (a.u.)'
%axis([0 (NN(1)-1)*dt*1000 0 max(abs([B1_init(:,1);B1_final(:,1)]))]);
axis([0 (NN(1)-1)*dt*1000 0 max(abs([B1_init(:,1);]))]);
subplot(222)
plot(0:dt*1000:(NN(1)-1)*dt*1000,angle(B1_init(1:NN(1),1)));
hold on
%plot(0:dt*1000:(NN(1)-1)*dt*1000,angle(B1_final(1:NN(1),1)),'r');
xlabel 'Time (ms)'
ylabel '\angle b_1(t) (Radians)'
axis([0 (NN(1)-1)*dt*1000 -pi pi]);
legend('Small-tip','Fast large-tip');

% meshplot Mz
subplot(223)
mesh(xx,yy,mzs0)
title 'Small-tip'
axis([-FOV/2 FOV/2 -FOV/2 FOV/2 -1 1]);
xlabel 'x (cm)'
ylabel 'y (cm)'
zlabel 'M_z'
subplot(224)
%mesh(xx,yy,mzs)
title 'Fast large-tip'
axis([-FOV/2 FOV/2 -FOV/2 FOV/2 -1 1]);
xlabel 'x (cm)'
ylabel 'y (cm)'
zlabel 'M_z'


% plot mag and phase, compare with target phase pattern
figure;
%subplot(1,2,1); imagesc(abs(m)-sin(tipangle/180*pi),[-0.1 0.1]);
subplot(1,2,1); imagesc(abs(m),[0 sin(tipangle/180*pi)*1.2]);
subplot(1,2,2); imagesc(angle(m./d),[-pi/4 pi/4]);

