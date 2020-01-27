clear; 
mse_blochAll = []; 
mse_projAll = []; 
pulse_durAll = [];
nrun = 3; 
for irun = 1:nrun
   % main file for continuous trajectory based on kt points
   close all;

   setupSimu_Cont;
   tipangle = 10;
   % slice_thickness = 0.5; %cm
   n_pe = 20;
   ishard = 1;
   Tread = 2.5e-3;
   
   % redo = 1; % 0: do not reselect kt points and connect them
   extype = 3; % 1: uniform excitation; 2: pre-phaseing; 3: cube
   if extype == 3
      d = zeros(size(b0map_simu));
      d(25:40,25:40,ceil(size(b0map_simu,3)/2)-3:ceil(size(b0map_simu,3)/2)+3) = 1;
      %     d(25:40,25:40,ceil(size(b0map_simu,3)/2)-4:ceil(size(b0map_simu,3)/2)) = 1;
      % d(25:40,25:40,:) = 1;
   elseif extype == 1
      d = ones(size(b0map_simu));
   elseif extype == 2
      d = exp(1i*b0map_simu*Tread*2*pi);
   end
   b0map = 1*b0map;
   
   [gbx gby gbz] = gradient(b0map);
   gb = gbx.^2+gby.^2 +gbz.^2;
   gb = ones(size(gb));
   %gb = smoothim(gb, 5, 2);
   %gb(:,:,1:6) = 0;
   %gb(:,:,11:end) = 0;
   weightIm = (gb).^0.5; % same dimension as image
   % weight = sparse(diag(weightIm(roi)));
   weight = diag_sp(weightIm(roi));
   
   save weighting weight weightIm gb;
   
   
   d_flip = d.*sin(tipangle./180*pi);
   d_flip(~mask) = 0;
   figure,im(d_flip); title('desired pattern (abs)');
   % [b1,gx,gy,gz] = spokerf(d,b0map,mask,tipangle,n_pe,slice_thickness,ishard,simuRange);
   % [b1,gx,gy,gz] = ktpointsRF(d,b0map,mask,tipangle,n_pe,ishard,simuRange);
   % [b1,gx,gy,gz] = ktContRF(d,b0map,mask,tipangle,n_pe,ishard,xyzRange);
   tolArray = [0:1:12];
   for it = 1:length(tolArray)
      it
      if it == 1
         redo = 1;
      else
         redo = 0;
      end
      tol = tolArray(it);
      [b1,gx,gy,gz,projIm,M_bloch] = ktContRF_greedy2(d,b0map,mask,tipangle,n_pe,ishard,simuRange,tol,redo);
      projImAll{it} = projIm;
      M_blochAll{it} = M_bloch;
      pulse_dur(it) = length(b1)*dt;
      figure,plot3(gx,gy,gz); 
   end
   
   datetime = datestr(now);
   dataname = [sprintf('errorVSsmooth_curve_%d_npe_%d',extype,n_pe) datetime];
   dataname(end-8) = [];
   dataname(end-5) = [];
   dataname(end-2) = [];
%    save(dataname);
   
   %% plot
   close all;
   % dataname = sprintf('errorVSsmooth_curve_%d',extype);
%    load(dataname);



   % load errorVSsmooth_curve.mat
   % projIm = projImAll{1};
   % M_bloch = M_blochAll{1};
   % diff_bloch = embed((M_bloch(mask)-d_flip(mask)),mask);
   % diff_proj = embed((projIm(mask)-d_flip(mask)),mask);
   % % diff_bloch = embed((abs(M_bloch(mask))-abs(d_flip(mask))),mask);
   % % diff_proj = embed((abs(projIm(mask))-abs(d_flip(mask))),mask);
   % % figure, im(d_flip,'cbar'); title('desired pattern');
   % figure, im(M_bloch,'cbar'); colormap default; title('bloch simu');
   % figure,im(diff_bloch,'cbar'); colormap default; title('excitation error');
   % figure,im(diff_proj,'cbar'); colormap default; title('excitation error (projection)');
   
   maxMask = mask_FOV(d);
   for it = 1:length(tolArray)
      
      projIm = projImAll{it};
      M_bloch = M_blochAll{it};
      diff_bloch = embed((M_bloch(mask)-d_flip(mask)),mask);
      if extype ~=3
      figure,im(diff_bloch,[0 max(abs(d_flip(:)))],'cbar'),colormap default;
      else
      figure,im(M_bloch,[0 max(abs(d_flip(:)))],'cbar'),colormap default; 
      end
      %     mse_proj(it) = norm(abs(projIm(mask))-abs(d_flip(mask)));
      %     mse_bloch(it) = norm(abs(M_bloch(mask))-abs(d_flip(mask)));
      %     max_proj(it) = max(abs(projIm(mask))-abs(d_flip(mask)));
      %     max_bloch(it) = max(abs(M_bloch(mask))-abs(d_flip(mask)));
      
      mse_proj(it) = norm(projIm(mask)-d_flip(mask))./norm(d_flip(mask));
      mse_bloch(it) = norm(M_bloch(mask)-d_flip(mask))./norm(d_flip(mask));
      max_proj(it) = max(abs(projIm(mask)-d_flip(mask)));
      max_bloch(it) = max(abs(M_bloch(mask)-d_flip(mask)));
      
      if extype == 3
         max_proj_out(it) = max(abs(projIm(~maxMask))-d_flip(~maxMask));
         max_bloch_out(it) = max(abs(M_bloch(~maxMask))-d_flip(~maxMask));
      end
   end
   % figure,plot(tolArray,mse_proj), hold on, plot(tolArray,mse_bloch,'r--'); xlabel('smoothing para (tol)'); ylabel('MSE');
   % figure,plot(pulse_dur,mse_proj), hold on, plot(pulse_dur,mse_bloch,'r--'); xlabel('pulse duration (ms)'); ylabel('MSE');
   % figure,plot(pulse_dur,max_proj), hold on, plot(pulse_dur,max_bloch,'r--'); xlabel('pulse duration (ms)'); ylabel('Max abs error');
   
   mse_blochAll = [mse_blochAll; mse_bloch];
   mse_projAll = [mse_projAll; mse_proj];
   pulse_durAll = [pulse_durAll; pulse_dur];
   diff_blochAll{it} = diff_bloch; 
   if extype == 3
      % plot the diff image in the out-ROI region
      diff_bloch = embed((abs(M_bloch(~maxMask))-d_flip(~maxMask)),~maxMask);
      figure,im(diff_bloch);
      figure,plot(pulse_dur,max_proj_out), hold on, plot(pulse_dur,max_bloch_out,'r--');
      xlabel('pulse duration (ms)'); ylabel('Max abs error');
   end
   
end
save(dataname); 
figure,plot(tolArray,mse_projAll,'b'), hold on, plot(tolArray,mse_blochAll,'r--'); xlabel('smoothing para (tol)'); ylabel('MSE');
figure,
for irun = 1:nrun
plot(pulse_durAll(irun,:),mse_projAll(irun,:),'b'), 
hold on, plot(pulse_durAll(irun,:),mse_blochAll(irun,:),'r--'); xlabel('pulse duration (ms)'); ylabel('MSE');
% hold on; plot(pulse_durAll,mse_projAll,'b'); hold on; plot(pulse_dur,mse_blochAll,'r--');
legend('proj','bloch'); 
end
figure,
for irun = 1:nrun
    hold on, plot(pulse_durAll(irun,:),mse_projAll(irun,:),'k'),
    xlabel('pulse duration (ms)'); ylabel('MSE');
    % hold on; plot(pulse_durAll,mse_projAll,'b'); hold on; plot(pulse_dur,mse_blochAll,'r--');
    % legend('proj','bloch');
end