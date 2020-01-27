function writewav_hao(wavname, b1, g, nom_fa)
% writewav_hao(wavname, b1, g, nom_fa)
% works good for stfr5.e
gx = g(:,1); 
gy = g(:,2); 
gz = g(:,3); 

[b1, gx, gy, gz] = sub_prepare_for_wavfile(b1,gx,gy,gz); %,0,0);
nom_bw = 2000; % bandwidth of tip-down sinc pulse (Hz)
gssamp = max(gz);
[paramsint16 paramsfloat] = myrfstat2(abs(b1(:,1)), nom_fa, gssamp, nom_bw);
writewav_1(wavname, sprintf('waves'), b1, [gx, gy, gz], paramsint16, paramsfloat); 

function [b1,gx,gy,gz] = sub_prepare_for_wavfile(b1,gx,gy,gz) %,tipupdelx,tipupdely)

% b1 = b1*1e4;     % convert to Gauss

% apply RF delay
%[b1,gx,gy,gz] = sub_rfdelay(b1,gx,gy,gz,tipupdelx,tipupdely);

% make smooth RF ramp
ramp = abs(b1(1))*[linspace(0,1,10)]';
% ramp = [-ramp/2; flipud(-ramp/2); ramp];

b1 = [ramp; b1];
gx = [0*ramp; gx];
gy = [0*ramp; gy];
gz = [0*ramp; gz];

% make length even
b1 = makeevenlength(b1);   % units: Gauss
gx = makeevenlength(gx);   % units: G/cm
gy = makeevenlength(gy);
gz = makeevenlength(gz);

% fix gradient directions
gx = -gx;  % confirmed sep9,2012
gy = gy;  % confirmed sep9,2012, update: on sep22,2013, seems should be negate as well 
gz = gz; % negate on sep22,2013, main_ktCont.m simulation is reversed to acquired image in z direction
%gx1tmp = gx1;
%gx1 = gy1; gy1 = gx1tmp; b11 = i*b11; %b11 = conj(b11);

% Conjugate to compensate for modulator phase conjugation. Confirmed Feb 2011.
% headcoil: +phase on theta channel produced -phase in image -> also conjugate for headcoil
b1 = conj(b1);

return;

function [paramsint16, paramsfloat] = myrfstat2(b1, nom_fa, g, nom_bw)
% Calculate RF parameters needed for RFPULSE struct in .e file.
% Needed for B1 scaling, SAR calculations, and enforcing duty cycle limits.
% See also mat2signa_krishna.m
%
% b1         real 1D vector containing B1 amplitude, size Nx1 [Gauss]
% nom_fa     nominal flip angle (degrees)
% g          slice-select gradient amp, G/cm (used for 2D multislice)
% nom_bw     Hz
% npre       delay between start of gradient waveform and start of b1 waveform


dt = 4e-6;                        % use 4 us RF sample width
%gamma = 4.2575e3;                  % Hz/Gauss
tbwdummy = 2;

% figure; plot(abs(b1));

% trim zeros before and after pulse (gradient rewinders are played here)
I = find(abs(b1)>1e-10);
b1 = b1(I(1):I(end));

npre = I(1)-1;
ndur = length(b1);
if mod(npre,2)
	npre = npre-1;
	ndur = ndur+1;
end
if mod(ndur,2)
	ndur = ndur+1;
end

fprintf(1,'npre, ndur = %d, %d \n', npre, ndur);

hardpulse = max(abs(b1)) * ones(length(b1),1);    % hard pulse of equal duration and amplitude

pw            = length(b1)*dt*1e3;                                       % ms
abswidth      = sum(abs(b1)) / sum(abs(hardpulse));
effwidth      = sum(b1.^2)   / sum(hardpulse.^2);
% or equivalently:  effwidth = sum(b1.^2)/(max(abs(b1))^2)/length(b1)
area          = abs(sum(b1)) / abs(sum(hardpulse)); 
dtycyc        = length(find(abs(b1)>0.2236*max(abs(b1)))) / length(b1);
maxpw         = dtycyc;
num           = 1;
%max_b1        = max(abs(b1)) ;                                            % Gauss
max_b1        = 0.2 ;                                                     % Gauss
max_int_b1_sq = max( cumsum(b1.^2)*dt*1e3 );                              % Gauss^2 - ms
max_rms_b1    = sqrt(mean(b1.^2));                                        % Gauss
nom_fa        = nom_fa;                                                   % degrees
nom_pw        = length(b1)*dt*1e6;                                        % us
%nom_bw        = tbwdummy / (dt * length(b1));                             % Hz
% max_int_b1    = abs(sum(b1))*dt*1000

% calculate equivalent number of standard pulses
stdpw = 1;                                                % duration of standard pulse (ms)
stdpulse = 0.117 * ones(round(stdpw/(dt*1e3)),1);
numstdpulses = num * effwidth * (pw/stdpw) * (max(abs(b1))/0.117)^2;


pulse50 = 0.117/180*50 * ones(round(stdpw/(dt*1e3)),1);
num50 = sum(abs(b1).^2)/sum(pulse50.^2)

paramsint16 = [npre ndur];
paramsfloat = [pw     abswidth  effwidth area          dtycyc      ...
               maxpw  num       max_b1   max_int_b1_sq max_rms_b1  ...
               nom_fa nom_pw    nom_bw   g numstdpulses              ];
               
return;

function writewav_1(fname,desc,B,g,paramsint16,paramsfloat)
% function writewav(fname,desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat)
% wav writting function for swstfr.e
% Single-shot rf waveforms for parallel transmission
%
% INPUT:
%   desc           ASCII description (arbitrary content and length).
%   rho            (Gauss) 
%                  size(rho) = N x Ncoils, where N = # waveform samples, 
%                  and Ncoils = number of transmit channels
%   theta          (radians). size(theta) = size(rho)
%   gx,gy,gz       (Gauss/cm). vector of length N
%   gmax           (Gauss/cm) maximum gradient amplitude supported by scanner  [4.0 G/cm]
%   paramsint16    Vector containing various RF pulse parameters needed for SAR and B1 scaling
%                  calculations -- just a placeholder for now. Max length = 32.
%   paramsfloat    Vector containing various RF pulse parameters needed for SAR and B1 scaling
%                  calculations -- just a placeholder for now. Max length = 32.
%
% Jon-Fredrik Nielsen, jfnielse@umich.edu
% $Id: writewav_hao.m,v 1.1 2015/05/22 21:39:26 jfnielse Exp $

rho = abs(B);
theta = angle(B);
gx = g(:,1);
gy = g(:,2);
gz = g(:,3);

% max length of params* vectors
nparamsint16 = 32;
nparamsfloat = 32;

max_b1 = 0.2;                 % fixed value, waveform is scaled relative to this

gmax  = 4.0;                  % Gauss/cm
% b1max = max(abs(rho(:)));     % Gauss

if (numel(paramsint16)>nparamsint16)
  error('writewav:nparamsint16',   'too many int16 parameters');
end
if (numel(paramsfloat)>nparamsfloat)
  error('writewav:nparamsfloat', 'too many float parameters');
end

% are waveforms the right size?
rhosize   = size(rho);
thetasize = size(theta);
L = [max(rhosize) max(thetasize) numel(gx) numel(gy) numel(gz)];
if any(L-max(rhosize)) 
  error('writewav:unequallength', 'all waveforms must have the same length');
end
if any(rhosize ~= thetasize) 
  error('writewav:unequallength', 'rf and theta waveforms must have the same size and dimensions');
end

% convert waveforms to column-major form
[nr,nc] = size(rho);    if nc > nr;   rho   = rho';     end;
[nr,nc] = size(theta);  if nc > nr;   theta = theta';   end;
[nr,nc] = size(gx);     if nc > nr;   gx    = gx';      end;
[nr,nc] = size(gy);     if nc > nr;   gy    = gy';      end;
[nr,nc] = size(gz);     if nc > nr;   gz    = gz';      end;

% remember to make sure waveforms have even length

% convert params* to row vectors, and pad to max length
paramsint16  = reshape(paramsint16,   1, numel(paramsint16));
paramsfloat  = reshape(paramsfloat, 1, numel(paramsfloat));
paramsint16  = [paramsint16   zeros(1, nparamsint16-numel(paramsint16))];
paramsfloat  = [paramsfloat zeros(1, nparamsfloat-numel(paramsfloat))];

[res,ncoils] = size(rho);


%
% write to file
%
fid = fopen(fname, 'w', 'ieee-be');

% write header
globaldesc = sprintf('RF waveform file for ssfpbanding project.\n');  
globaldesc = sprintf('%sCreated by %s.m on %s.\n', globaldesc, mfilename('fullpath'), datestr(now));  
globaldesc = sprintf('%snchannels = %d, res = %d\n', globaldesc, ncoils, res);  
desc = sprintf('%s%s\n', globaldesc, desc);
fwrite(fid, numel(desc), 'int16');      % number of characters in ASCII description
fwrite(fid, desc, 'uchar');

fwrite(fid, ncoils, 'int16');          % shorts must be written in binary -- otherwise it won't work on scanner 
fwrite(fid, res,    'int16');
fprintf(fid, 'b1max:  %f\n', max_b1);   % floats are OK in ASCII on scanner 
fprintf(fid, 'gmax:   %f\n', gmax);

fwrite(fid, nparamsint16, 'int16');
fwrite(fid, paramsint16,  'int16');
fwrite(fid, nparamsfloat, 'int16');
for n = 1:nparamsfloat
	fprintf(fid, '%f\n', paramsfloat(n)); 
end

% write binary waveforms (*even* short integers -- the psd sets the EOS bit, so don't have to worry about it here)
maxiamp = 2^15-2;                           % max instruction amplitude (max value of signed short)
rho   = 2*round(rho/max_b1*maxiamp/2);
theta = 2*round(theta/pi*maxiamp/2);
gx    = 2*round(gx/gmax*maxiamp/2);
gy    = 2*round(gy/gmax*maxiamp/2);
gz    = 2*round(gz/gmax*maxiamp/2);
fwrite(fid, rho(:),   'int16');
fwrite(fid, theta(:), 'int16');
fwrite(fid, gx(:),    'int16');
fwrite(fid, gy(:),    'int16');
fwrite(fid, gz(:),    'int16');

fclose(fid);

% EOF







