% main_opt
close all; 
setupSimu_Cont; 
tipangle = 10; 
ishard = 1;
Tread = 2.5e-3;
extype = 3; % 1: uniform excitation; 2: pre-phaseing; 3: rectangle 
if extype == 3
    n_pe = 40; 
    d = zeros(size(b0map_simu));
    d(25:40,25:40,ceil(size(b0map_simu,3)/2)-3:ceil(size(b0map_simu,3)/2)+3) = 1;
%     d(25:40,25:40,ceil(size(b0map_simu,3)/2)-4:ceil(size(b0map_simu,3)/2)) = 1;
    % d(25:40,25:40,:) = 1;
elseif extype == 1
    d = ones(size(b0map_simu));
    d(~mask) = 0; 
elseif extype == 2
    n_pe = 5; 
    d = exp(1i*b0map_simu*Tread*2*pi); 
    d(~mask) = 0; 
end
d_flip = sin(tipangle/180*pi)*d; 
Trf = 2e-3; 
g0 = zeros(15,1); 

% [X,FVAL] = fmincon(@(g) func_of_g(g,d,b0map,roi,sens,xyzRange,Trf), g0, [],[],[],[],-1*ones(size(g0)),1*ones(size(g0))); 
% save opt_all
% 
% % the following is only avialable on matlab 2013
% opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter','PlotFcns',@optimplotfval);
% problem = createOptimProblem('fmincon','objective',...
%  @(g) func_of_g(g,d_flip,b0map,roi,sens,xyzRange,Trf),'x0',g0,'lb',-4*ones(size(g0)),'ub',4*ones(size(g0)),'options',opts);
% gs = GlobalSearch; 
% [X,FVAL,EXITFLAG,OUTPUT] = run(gs,problem)

niter = 100; 
options = optimset('Algorithm','interior-point','MaxIter',niter); % run interior-point algorithm
[X,FVAL,EXITFLAG,OUTPUT] = ...
   fmincon(@(g) func_of_g(g,d_flip,b0map,roi,sens,xyzRange,Trf), g0, [],[],[],[],-4*ones(size(g0)), 4*ones(size(g0)), [], options); 
save opt_all_global

load opt_all_global
% X = zeros(size(X)); 
[b1,gx,gy,gz,projIm]= func_of_g_blochIm(X, d_flip,b0map,roi,sens,xyzRange,Trf); 
T1 = 1000; T2 = 100; dt = 4e-6;
M_pr = parallel_blochCim(0,b1,gx,gy,gz,sens,xyzRange.x,xyzRange.y,xyzRange.z,dt,b0map,roi,T1/1000,T2/1000);
M_diff = M_pr - d_flip; 
M_diff_pro = projIm - d_flip;
nrmse = norm(M_diff(:))/norm(d_flip(:)); 
figure,im(abs(M_diff_pro),[0 sin(tipangle/180*pi)],'cbar'),colormap default; 
figure,im(abs(M_diff),[0 sin(tipangle/180*pi)],'cbar'),colormap default; 
title(sprintf('nrmse: %f',nrmse)); 
figure,im(angle(M_pr./d)*180/pi,'cbar'),colormap default; 
