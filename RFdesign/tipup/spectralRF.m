function [b1,mxy] = spectralRF(Trf,TE,wn,signOfTargetPhase,lambda,type,varargin)
% function spectralRF(Trf,TE,signOfTargetPhase,lambda,type,varargin)
%
% Inputs
%   Trf                 [1 1]       pulse duration (sec)
%   TE                  [1 1]       target TE (sec)
%   wn                  [1 nfreq]   target frequency range (Hz)
%   signOfTargetPhase   +1 or -1    For minimal T2*-weighting (bSSFP-like contrast) when used as tip-up pulse in STFR, set to +1
%   lambda              [1 1]       Tikhonov regularization constant (0.6 or 1.0 seems to work well)
%   type                [string]    'tipdown' or 'tipup'
% Options:
%   Tfree     [1 1]     free precession time (sec). Default: Tfree = TE.
%   fmt       [string]  plot formatting string. Default: 'b'.

% parse input options
arg.Tfree = TE;
arg.fmt = 'b';
arg.hold = 'off';
arg = toppe.utils.vararg_pair(arg, varargin);
Tfree = arg.Tfree;
fmt = arg.fmt;

if Trf > 20e-3
	error('Trf is in sec');
end
if TE > 50e-3
	error('TE is in sec');
end

if strcmp(type, 'tipdown')
	signOfB0DuringRF = 1;
else
	signOfB0DuringRF = -1;   % because we're playing the tip-up pulse in time-reverse
end

if abs(signOfTargetPhase) ~= 1
%	error('signOfTargetPhase must be 1 or -1');
end

%% Target frequency range
%wn = [-20:0.1:20]';             % Hz
wn = wn(:);

%% Target frequency response
flip = 15;      % degrees
d = sind(flip)*exp(signOfTargetPhase*1i*2*pi*wn*TE); 
       
%% Construct A
m0 = 1; 
gambar = 4257; %Hz/G
gam = gambar*2*pi;
dt = 4e-6; %s; 
t = dt:dt:Trf;
t = t-Trf;
[b0tt, tt] = ndgrid(signOfB0DuringRF*wn(:),t(:));
A = 1i*gam*m0*exp(1i*2*pi*(b0tt.*tt))*dt;

%% Solve for b1   
%lambda = 0.6; % 0.6 for nov11 jon
b1 = pinv(A'*A+lambda^2*eye(length(t),length(t)))*A'*d; 

% % penalize the pulse more at the beginning to reduce the over tipping in the process   
% w = linspace(2,0.7,length(t)); W = diag(w);  
% b1 = pinv(A'*A+W.*W*eye(length(t),length(t)))*A'*Mxy; 

%% Simulate
T1 = 1000; T2 = 80;    % msec
dt = 4e-3;             % msec
w = 1*wn(:);           % simulate over this frequency range

nfreq = length(w);

b1t = b1(:)*1e-4;       % [nstep 1], Tesla
nstep = length(b1t);

if strcmp(type, 'tipup');
	% magnetization at end of Tfree
	m0 = [sind(flip)*exp(1i*2*pi*w*Tfree) 0*ones(size(w)) cosd(flip)*ones(size(w))]; 

	% scale b1 to match the (on-resonance) tipdown flip angle
   Beff = [real(b1t) imag(b1t) 0*ones(nstep,1)];     % [nstep 3], Tesla
   mtmp = toppe.utils.rf.blochsim([0 0 1], Beff, T1, T2, dt, nstep);
   mxy = mtmp(end,1) + 1i*mtmp(end,2);
	b1t = b1t*sind(flip)/abs(mxy);
	fprintf('pulse area/sar = %.2f/%.2f x standard pulse \n', sum(abs(b1t*1e4))/(250*0.117), sum(abs(b1t*1e4).^2)/sum(.117^2*ones(250,1)));

	% time-reverse and negate
	b1t = flipud(-b1t);   
else
	fprintf('pulse area/sar = %.2f/%.2f x standard pulse \n', sum(abs(b1t*1e4))/(250*0.117), sum(abs(b1t*1e4).^2)/sum(.117^2*ones(250,1)));
	% start with longitudinal magnetization
	m0 = [0*ones(size(w)) 0*ones(size(w)) ones(size(w))];
end

% effective field (z component)
Bz = ones(nstep,1)*w(:)'/gambar*1e-4;     % [nstep nfreq], Tesla

for ii = 1:nfreq
   Beff = [real(b1t) imag(b1t) Bz(:,ii)];     % [nstep 3], Tesla
   mtmp = toppe.utils.rf.blochsim(m0(ii,:), Beff, T1, T2, dt, nstep);
   mxy(ii) = mtmp(end,1) + 1i*mtmp(end,2);
	mz(ii) = mtmp(end,3);
end

%% Display
T = dt*(1:nstep);
subplot(151); hold(arg.hold); plot(T,abs(b1),fmt); xlabel('time (ms)'); ylabel('abs(b1), Gauss'); title('abs(b1)');

subplot(152); hold(arg.hold); plot(T,angle(b1),fmt); title('angle(b1)');

subplot(153); hold(arg.hold); plot(w,abs(mxy),fmt); hold on; plot(w, sind(flip)*ones(size(w)), 'r--'); 
xlabel('frequency (Hz)'); ylabel('abs(mxy)'); axis([w(1) w(end) 0.0 1.0]);
title(['Mxy at end of ' type ' pulse']);

subplot(154); hold(arg.hold); plot(w,angle(mxy)/pi*180,fmt); hold on; plot(w, signOfB0DuringRF*signOfTargetPhase*w*TE*180, 'r--');
xlabel('frequency (Hz)'); ylabel('angle(mxy), degrees'); axis([w(1) w(end) -100 100]);
title(['angle(Mxy) at end of ' type ' pulse']);

subplot(155); hold(arg.hold); plot(w,mz,fmt); xlabel('frequency (Hz)'); ylabel('mz'); axis([w(1) w(end) -1 1]);
title(['Mz at end of ' type ' pulse']);


