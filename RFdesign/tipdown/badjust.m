%function M_bloch = badjust
% Do additive angle correction, to suppress OV excitation

fprintf(1,'Remember to load jfntmp (created in ktCont_greedy2.m), b0map_simu, and d_flip \n');
%load jfntmp  % see ktCont_greedy2.m
load b0map_simu
load d_flip;

ov_roi2 = ov_roi;
ov_roi2(:,:,10:16) = 0;
weightIm2 = weightIm;
weightIm2(~ov_roi2) = 0.03;
W2 = diag_sp(weightIm2(mask));

% form A and get b1
dt = 4e-6; tt = [0:length(K_rv)-1]*dt-length(K_rv)*dt;
As = formA([K_rv(:,1)*xfov, K_rv(:,2)*yfov, K_rv(:,3)*zfov], permute(sens,[2,3,4,1]), mask, fieldmap, tt, xyzRange);
betavec = sqrt(beta)*ones(size(As,2),1);
%weightImaa = weightIm;
%weightImaa(ov_roi) = 5*weightImaa(ov_roi); % add extra weight tOV spins
W = diag_sp(weightIm(mask));
%b1 = bupdate(zeros(size(betavec)), As, W, d, mask, betavec, 100);
%b1 = b1.*sin(tipangle/180*pi);

% get small-tip exictation (projIm)
g = k2g_hao(K_rv);
gx = g(:,1);
gy = g(:,2);
gz = g(:,3);
projIm = embed(As*b1,mask);

% Compare with Bloch simulation
T1 = 1000; T2 = 100; dt = 4e-6;
b0map_simu_zero = 0*b0map_simu;
Nt = length(gx);
M_pr = parallel_blochCim(0,double(reshape(b1,[Nt,ncoils])),double(gx),double(gy),double(gz),...
                        double(sens),double(xyzRange.x),double(xyzRange.y),double(xyzRange.z),...
                        dt,double(b0map_simu),mask,T1/1000,T2/1000);

% projIm and M_pr differ, presumably due to violations of small-tip assumption.
% Now we try to reduce OV excitation using Grissom's additive angle method for large-tip RF pulse design.
load rois;   % need ov_roi
mask_ov = mask;
mask_ov(~ov_roi) = 0;
W_ov = diag_sp(weightIm(mask_ov));
A_ov = formA([K_rv(:,1)*xfov, K_rv(:,2)*yfov, K_rv(:,3)*zfov], permute(sens,[2,3,4,1]), mask_ov, fieldmap, tt, xyzRange);
betavec_ov = sqrt(100*beta)*ones(size(A_ov,2),1);
mv_in = M_pr(ov_roi2 & mask);
Nov = numel(mv_in);
nrmse_ov_in = norm(mv_in)/sqrt(Nov)/max(abs(d_flip(:)));
M_pr2 = M_pr;
b1new = b1;
betavec = sqrt(100)*ones(size(As,2),1);
clear nrmse_ov_out;
for iaa=1:3
	deltad = -M_pr2; %/max(abs(M_pr2(:)));  % target pattern for deltab1 (correction pulse)
	deltad(~ov_roi2) = 0;
	%deltab1 = bupdate(zeros(size(betavec_ov)), A_ov, W_ov, deltad, mask_ov, betavec_ov, 40); 
	%deltab1 = bupdate(deltab1, A_ov, W_ov, deltad, mask_ov, betavec_ov, 10); 
	deltab1 = bupdate(zeros(size(betavec)), As, W2, deltad, mask, betavec, 30); 
	b1new = b1new+deltab1;
	M_pr2 = parallel_blochCim(0,double(reshape(b1new,[Nt,ncoils])),double(gx),double(gy),double(gz),double(sens),double(xyzRange.x),double(xyzRange.y),double(xyzRange.z),dt,double(b0map_simu_zero),mask,T1/1000,T2/1000);
	mv_out = M_pr2(ov_roi2 & mask);
	nrmse_ov_out(iaa) = norm(mv_out)/sqrt(Nov)/max(abs(d_flip(:)))
end

mv_in = M_pr(ov_roi2 & mask);
mv_out = M_pr2(ov_roi2 & mask);
Nov = numel(mv_in);
nrmse_ov_in = norm(mv_in)/sqrt(Nov)/max(abs(d_flip(:)))

return;

% try direct excitation design -- no doesn't work
%{
d2 = d_flip;
d2(ov_roi) = -M_pr(ov_roi);
b1new2 = bupdate(b1, As, W, d2, mask, betavec, 50); 
M_pr2 = parallel_blochCim(0,double(reshape(b1new2,[Nt,ncoils])),double(gx),double(gy),double(gz),double(sens),double(xyzRange.x),double(xyzRange.y),double(xyzRange.z),dt,double(b0map_simu_zero),mask,T1/1000,T2/1000);
%}

figure; plot(nrmse_ov_out);
figure; im(cat(1,ov_roi.*projIm,ov_roi.*M_pr,ov_roi.*M_pr2)); colormap default;
figure; im(cat(1,abs(projIm-d_flip),abs(M_pr-d_flip),abs(M_pr2-d_flip)),[0 0.1]); colormap default;
figure; 
subplot(121); hold on; z=7; y=32; plot(abs(M_pr(:,y,z)),'r'); plot(abs(M_pr2(:,y,z)),'b'); xlabel('distance (pixels)'); ylabel('magnetization (a.u.)'); legend('before correction','after additive-angle correction');
subplot(122); hold on; z=7; x=32; plot(squeeze(abs(M_pr(x,:,z))),'r'); plot(squeeze(abs(M_pr2(x,:,z))),'b'); xlabel('distance (pixels)'); ylabel('magnetization (a.u.)'); legend('before correction','after additive-angle correction');

% write pulses before and after correction to .wav files
mat2wav(abs(b1),angle(b1),gx,gy,gz,tipangle,'badjust_before.wav','IV pulse');
mat2wav(abs(b1new),angle(b1new),gx,gy,gz,tipangle,'badjust_after.wav','IV pulse');


