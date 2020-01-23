function [b1,mxy] = spectralRF_jfn(Trf,TE,signOfTargetPhase,lambda,type,varargin)
% function spectralRF_jfn(Trf,TE,signOfTargetPhase,lambda,type,varargin)
%
% Inputs
%   Trf                 [1 1]     pulse duration (sec)
%   TE                  [1 1]     target TE (sec)
%   signOfTargetPhase   +1 or -1
%   lambda              [1 1]     Tikhonov regularization constant (0.6 or 1.0 seems to work well)
%   signOfB0DuringRF    [1 1]   +1 or -1
%   type                [string]  'tipdown' or 'tipup'
% Options:
%   Tfree     [1 1]     free precession time (sec). Default: Tfree = TE.
%   fmt       [string]  plot formatting string. Default: 'b'.

% Results/notes:
%  See also Google doc 'Figures for weekly meeting w/ Doug' (jfnielse@umich.edu)
%
%  >> spectralRF_jfn(6e-3,5e-3, sgnTgtPhs=1,lam=1.0,sgnB0DurRF=1);  % phase slope exceeds target phase; flip below target near 0 Hz, increases rapidly as |w| increases
%  >> spectralRF_jfn(5e-3,5e-3,-1,1.0,-1);    % phase slope exceeds target phase; flip below target near 0 Hz, increases rapidly as |w| increases
%
% Compare
%  >>  spectralRF_jfn(2e-3,4e-3,-1,0.999,'tipup',6e-3);   % this produces greater signal suppression away from w=0 => should enhance BOLD contrast
%  >>  spectralRF_jfn(2e-3,0e-3,-1,0.999,'tipup',6e-3);
% 
% Potentially useful designs:
%  >> spectralRF_jfn(2e-3,5e-3,1,1.4,'tipdown');  % matches target phase over [-20:20] Hz reasonably well
%  >> spectralRF_jfn(3e-3,5e-3,-1,0.999,'tipup',10e-3);   % creates a sharper freq profile (after tip-up) than a hard tipup
%  or
%  >> spectralRF_jfn(2e-3,5e-3,-1,1.2,'tipup',5e-3);  % pulse area/sar = 0.58/0.22 x standard pulse

%

% TE = 8e-3;
% Trf = 3e-3;

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
	error('signOfTargetPhase must be 1 or -1');
end

%% Target frequency range
wn = [-20:0.1:20]';             % Hz

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

%subplot(121); plot(abs(b1));
%subplot(122); plot(angle(b1));

% % penalize the pulse more at the beginning to reduce the over tipping in the process   
% w = linspace(2,0.7,length(t)); W = diag(w);  
% b1 = pinv(A'*A+W.*W*eye(length(t),length(t)))*A'*Mxy; 

% % penalize the pulse more at the beginning to reduce the over tipping in the process
% w = linspace(3*lambda,lambda,length(t)); W = diag(w);
% b1 = pinv(A'*A+W.*W*eye(length(t),length(t)))*A'*Mxy;

%% Simulate
T1 = 1000; T2 = 80;    % msec
dt = 4e-3;             % msec
w = wn(:);           % simulate over this frequency range

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
subplot(151); hold(arg.hold); plot(T,abs(b1),fmt); xlabel('time (ms)'); ylabel('abs(rf), Gauss');
subplot(152); hold(arg.hold); plot(T,angle(b1),fmt); 
subplot(153); hold(arg.hold); plot(w,abs(mxy),fmt); hold on; plot(w, sind(flip)*ones(size(w)), 'r--'); 
xlabel('frequency (Hz)'); ylabel('abs(mxy)'); axis([w(1) w(end) 0.0 1.0]);
subplot(154); hold(arg.hold); plot(w,angle(mxy)/pi*180,fmt); hold on; plot(w, signOfB0DuringRF*signOfTargetPhase*w*TE*180, 'r--');
xlabel('frequency (Hz)'); ylabel('angle(mxy), degrees'); axis([w(1) w(end) -100 100]);
subplot(155); hold(arg.hold); plot(w,mz,fmt); xlabel('frequency (Hz)'); ylabel('mz'); axis([w(1) w(end) -1 1]);


