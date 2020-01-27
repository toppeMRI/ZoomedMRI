% creat sinc wavfile
% sinc_main.m
% not work on the scanner. Not sure why? 

%  Hao Sun : Dec15,2012; the university of Michigan
clear
if exist('myslr') ~= 2
    dosetup; 
end
% doSimu = 0; 
% setupPara; % load maps and setup parameters
% load designParas; 

%% create tip-down pulse. slice-selective.

opslthick = 16;  % cm
nom_fa = 35; 
tbw = 4; 
tipdowndur = 1.5;        % msec
Tsinc = 1.11;       % time from middle of sinc to end of tip-down gz waveform, EXCLUDING gzdep
isbssfp = 1; 
ncycles = 12/1*opslthick; % crusher gradient: 12 cycle/cm
[b1,gz,gz2refocus] = myslrrf(opslthick,nom_fa,tbw,tipdowndur,'ls',isbssfp,ncycles);
gx = zeros(size(gz));
gy = zeros(size(gz));
if isbssfp == 0
   if opslthick < 1
      writewav_hao(sprintf('sinc%d_%dmm.wav',nom_fa, opslthick*10),b1,[gx,gy,gz],nom_fa);
   else
      writewav_hao(sprintf('sinc%d_%dcm.wav',nom_fa, opslthick),b1,[gx,gy,gz],nom_fa);
   end
else
   if opslthick < 1
      writewav_hao(sprintf('bssfp%d_%dmm.wav',nom_fa, opslthick*10),b1,[gx,gy,gz],nom_fa);
   else
      writewav_hao(sprintf('bssfp%d_%dcm.wav',nom_fa, opslthick),b1,[gx,gy,gz],nom_fa);
   end
end
   
return; 

% %% write to .wav file that can be loaded into the psd (stfr.e)
% [b1, gx, gy, gz] = sub_prepare_for_wavfile(b1,gx,gy,gz);
% nom_bw = 2000; % bandwidth of tip-down sinc pulse (Hz)
% gssamp = max(gz); 
% [paramsint16 paramsfloat] = myrfstat2(abs(b1(:,1)), nom_fa, gssamp, nom_bw);
% wavfileName = ['sinc_' num2str(nom_fa) '_' num2str(opslthick) 'cm.wav']; 
% writewav_1(wavfileName,sprintf('bssfp pulse, sinc'),b1,[gx,gy,gz],paramsint16,paramsfloat);
% plotwav(wavfileName);
% max(abs(b1)); 
