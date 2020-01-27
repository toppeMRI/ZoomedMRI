function [rf,gss,gzrep,npre,nrf] = myslrrf(opslthick,nom_flip,tbw,dur,type,isbalanced,ncycles)
% function [rf,gss,gzrep] = myslrrf(opslthick,nom_flip,tbw,dur,type,isbalanced)
%
% copied from ~/projects/psd/spiral3drf/rfgen/ on 08/29/11
% Create slice-selective SLR pulse, including gz gradient.
% Calls ~/matlab/rf/myslr.m.
%
% INPUTS:
%  opslthick  - cm 
%  nom_flip   - degrees
%  tbw        - total # zero crossings = time-bandwidth product
%  dur        - duration (msec)
%  type       - 'min' for minimum-phase design, 'ls' otherwise
%
% OUTPUTS:
%  rf         - b1 waveform (real), in Gauss
%  gz         - slice-select gradient waveform (Gauss/cm)
%  gzrep      - slice-select gradient refocusing lobe
%
% $Id: myslrrf.m,v 1.1 2015/05/22 21:39:24 jfnielse Exp $
%
%addpath('~/matlab/rf');
%addpath('~/matlab/rf_tools');

% rf = myslr(256,4,30,10,'min','tmp');

if ~exist('type')
	type = 'ls';
end

dt = 4e-3;                   % sample size (msec)
npix = round(dur/dt);

rf = myslr(npix,tbw,nom_flip,opslthick,type,'tmp');    % row vector

% scale to Tesla, and make duration even
rf = nom_flip/180*pi * rf / sum(rf);    % normalized such that flip angle (radians) = sum(rf)
gamma = 4.2575;                         % kHz/Gauss
rf = rf / gamma / dt /2/pi;             % since rf = dth = 2 *pi * gamma * b1 (in Gauss) * dt
% rf = rf*1e-4;                           % Tesla
rf = [rf(:); zeros(mod(length(rf),2),1)];
npix = length(rf);

% make slice-select gradient waveform
bw = tbw / dur;                    % kHz
g = bw / (gamma * opslthick);      % slice-select gradient amplitude (G/cm)
mxg = 3.9;                         % G/cm
mxs = 14.9;                        % G/cm/msec
if g > mxg
	error('g > max g');
end

% slice-select trapezoid 
gss = g*ones(1,npix);                             % plateau of slice-select gradient
s = mxs * dt;                                     % max change in g per sample (G/cm)
ramp = fliplr([0:s:(g-s)]);
ramp = [ramp zeros(1,mod(length(ramp),2)) ];      % makes res_kpre even

% slice-select rephaser 
gssarea = sum(gss) * dt * 1e-3 ;                           % G/cm*s
ramparea = sum(ramp) * dt * 1e-3 ;
gzrep = -trapwave((gssarea/2+ramparea), dt*1e-3, mxg, mxs*1000);
gzrep = [gzrep zeros(1,mod(length(gzrep),2)) ];            % make length even

% spoiler gradient (or dephaser for bSSFP)
if isbalanced
	gzdep = [fliplr(gzrep) fliplr(ramp)];
else
	%area = sum([gss]) * dt * 1e-3 ;                        % G/cm*s
	gamma = 4257.6;                                        % Hz/Gauss
	areatot = ncycles/(gamma*opslthick);                   % Desired total spoiling area, G/cm*sec
	gzdep = mybridged(areatot-gssarea/2,0,gss(1));
	gzdep = gzdep(:)';
end

% put together the complete slice-select gradient waveform
gss = [gzdep gss ramp gzrep]';

npre = length(gzdep);
nrf = length(rf);

% make gss and rf the same length
rf = [zeros(length(gzdep),1); rf; zeros(length(ramp)+length(gzrep),1)];

% plot
subplot(1,2,1); plot(abs(rf));
legend('abs(rf) (Gauss)');
subplot(1,2,2); plot(gss);
legend('gradient (Gauss/cm)');

return;
