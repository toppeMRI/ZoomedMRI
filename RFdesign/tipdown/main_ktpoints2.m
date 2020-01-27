
% main file for kt-points
close all;
clear;
%setupSimu_ktpoints; 
setupSimu_Cont; 
tipangle = 10; 
% slice_thickness = 0.5; %cm
n_pe = 25; 
ishard = 1;
Tread = 2.5e-3;
extype = 3; % 1: uniform excitation; 2: pre-phaseing; 3: rectangle 
if extype == 3
    d = zeros(size(b0map_simu));
    %d(25:40,25:40,ceil(size(b0map_simu,3)/2)-3:ceil(size(b0map_simu,3)/2)+3) = 1;
    d(25:40,25:40,ceil(size(b0map_simu,3)/2)-1:ceil(size(b0map_simu,3)/2)+3) = 1; 

%     d(25:40,25:40,:) = 1;
elseif extype == 1
    d = ones(size(b0map_simu)); 
elseif extype == 2
    d = exp(1i*b0map_simu*Tread*2*pi); 
end
b0map = 1*b0map; 
% [b1,gx,gy,gz] = spokerf(d,b0map,mask,tipangle,n_pe,slice_thickness,ishard,simuRange); 
[b1,gx,gy,gz] = ktpointsRF(d,b0map,mask,tipangle,n_pe,ishard,simuRange); 
% [b1,gx,gy,gz] = ktContRF(d,fieldmap,mask,tipangle,n_pe,ishard,xyzRange); 
% [b1,gx,gy,gz] = ktContRF_greedy(d,b0map,mask,tipangle,n_pe,ishard,simuRange); 
writewav_hao('ktpoints.wav',b1,[gx,gy,gz],tipangle);


