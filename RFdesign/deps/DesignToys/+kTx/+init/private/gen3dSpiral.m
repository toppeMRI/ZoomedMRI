%==========================================================================
% 11-1-10: script for generating 3D spiral. This performs smooth
% interpolation during the lead-out phase and eddy current preemphasis

% modified by Hao: convert unit to our convention
% 2019-03-15
% Ray: modified from Hao's generate_3d_spiral.m
function [g,k,S,g_old] = gen3dSpiral(varargin)
%% constant definitions
ntrail=0;
doPlots = false;  % plot

%% get inputs
arg = envMR('get_s'); % {'dt','gam','gMax','sMax'} {Sec, Hz/G, G/cm, G/cm/Sec}
arg.xlead = false;  % xmit on lead
arg.precomp = false;  %  preemphasis

arg.dur = 2.2e-3; % core duration

% kspace pars
arg.kMax = 12*pi;
arg.thMax = 20*pi;
arg.phiMax = 2*pi;
arg.R = [10 0.5];

arg = attrParser(arg, varargin);
[xlead, precomp, dur, kMax, thMax, phiMax, R, sMax, dt, gam] = getattrs(arg ...
  ,{'xlead','precomp','dur','kMax','thMax','phiMax','R','sMax','dt','gam'});

gamma_mT = gam*2*pi*10; % Rad Hz/mT

%% parameters + starting k-space trajectory
M=fix(dur/dt);

r=linspace(1,0,M).^1; %spiral in

% R(1)=0 is the flag for a linear trajectory
if R(1)~=0
  r=1./(exp(R(1)*(R(2)-r))+1);
end
kr = kMax*r;
th=linspace(0,thMax,M);
ph=linspace(0,phiMax,M);

% --- make kr go nicely to zero
nfit = fix(M*0.2); % arbitrarily use 20% samples
nfix = ceil(M*0.02);
idx = [(M-nfit)+(1:nfix) (M-nfix+1):M];
kr((M-nfix+1):M)=0;
kr_end = spline(idx,kr(idx),(M-nfit+1):M);
kr((M-nfit+1):M)=kr_end;
%-------------------------

kx = kr.*sin(th).*sin(ph);
ky = kr.*sin(th).*cos(ph);
kz = kr.*cos(th);
%k.list = [kx.' ky.' kz.' (t-T).'];

%%  initial gradient trajectory
kG = -gamma_mT*dt*triu(ones(M));
G = kG\[kx.' ky.' kz.'];

%% smooth lead-out
% # of smaples to fit inside gradient
nfit = 8;
% pad leading samples with a few zeros (for PCG)
npad = 5;
% minimum # samples for lead out
t_lead = max(G(1,:))/sMax;
nlead = ceil(t_lead/dt);

% fit for spline given initial zeros and then actual demand waveform
nn=nfit+npad+nlead;
tfit = dt*(1:nn);
idx = [1:npad (npad+nlead+1):nn];
Glead = zeros([nlead+npad 3]);
for ii=1:3
  gg = [zeros([npad 1]);G(1:nfit,ii)];
  gp = spline(tfit(idx),gg,tfit);
  Glead(:,ii) = gp(1:(nlead+npad));
end

%% trailing samples
% default is zero
Gtrail = zeros(ntrail,3);
g.Ntrail=ntrail;

%%  put together
if ~xlead % leading samples are separate
  g.Nlead = size(Glead,1);
  g.Ncore = size(G,1);
  tmp = [G; Gtrail];
  g.x = tmp(:,1);g.y = tmp(:,2);g.z = tmp(:,3);
  tmp = [Glead; G; Gtrail];
  g.all.x = tmp(:,1);g.all.y = tmp(:,2);g.all.z = tmp(:,3);
else % put leading samples in core section
  g.Nlead = 0;
  G = cat(1,Glead,G);
  g.Ncore = size(G,1);
  tmp = [G; Gtrail];
  g.x = tmp(:,1);g.y = tmp(:,2);g.z = tmp(:,3);
  g.all.x = tmp(:,1);g.all.y = tmp(:,2);g.all.z = tmp(:,3);
end

%% pre-emphasis
g_old = g;
if precomp
  g = grad_precomp(g);
end

%% Slew rate (based on pre-emphasized waveforms)
S = cat(2,gradient(g.all.x,dt),gradient(g.all.y,dt),...
  gradient(g.all.z,dt));

%% get k-space (based on demand waveforms)
M = size(G,1);
kG = -gamma_mT*dt*triu(ones(M));
kk = kG*G;

tt = dt*(0:M-1);
T=tt(end);

% interpolate in case dt_design~=dt_system
ti = 0:dt:T;
ki = interp1(tt,kk,ti);

kw = ti-ti(end);

% struct for pulse design
k.list = cat(2,ki,kw.');

%% get k space based on actual waveforms
Gall = [g.all.x, g.all.y, g.all.z];
Mall = size(Gall,1);
kGall = -gamma_mT*dt*triu(ones(Mall));
kkall = kGall*Gall;

ttall = dt*(0:Mall-1);
Tall=ttall(end);

% interpolate in case dt_design~=dt_system
tiall = 0:dt:Tall;
kiall = interp1(ttall,kkall,tiall);

kwall = tiall-tiall(end);

% struct for pulse design
k.all = cat(2,kiall,kwall.');

%% convert unit from Malik's code to our convention
% convert k from rad/m to cycle/cm
k.list = k.list/(2*pi)/100; %conver from rad/m to cycle/cm
k.all = k.all/(2*pi)/100; %conver from rad/m to cycle/cm

% convert g from mt/m to Gauss/cm
g.x = g.x/10;
g.y = g.y/10;
g.z = g.z/10;
g.all.x = g.all.x/10;
g.all.y = g.all.y/10;
g.all.z = g.all.z/10;

g_old.x = g_old.x/10;
g_old.y = g_old.y/10;
g_old.z = g_old.z/10;
g_old.all.x = g_old.all.x/10;
g_old.all.y = g_old.all.y/10;
g_old.all.z = g_old.all.z/10;

% convert slew rate S from mt/m/s to Gauss/cm/s
S = S/10;


%% plot gradients
if doPlots
  figure;
  set(gcf,'position',[360   187   368   735])
  nr=3;nc=1;
  lw=1.5;
  
  subplot(nr,nc,1);
  plot3(k.list(:,1),k.list(:,2),k.list(:,3),'linewidth',lw);
  grid
  
  tt = dt*(0:length(g.all.x)-1)*1e3;
  subplot(nr,nc,2);
  plot(tt(1:g.Nlead),[g.all.x(1:g.Nlead) g.all.y(1:g.Nlead) g.all.z(1:g.Nlead)]...
    ,'linewidth',lw/2);
  hold
  plot(tt(g.Nlead+1:end),[g.all.x(g.Nlead+1:end) g.all.y(g.Nlead+1:end) ...
    g.all.z(g.Nlead+1:end)] ,'linewidth',lw);
  legend('gx','gy','gz');
  grid
  xlabel('t/ms');ylabel('G : g cm^{-1}')%ylabel('G / mT m^{-1}')
  xlim(tt([1 end]))
  
  subplot(nr,nc,3);
  plot(tt,S*1e-3,'linewidth',lw);
  grid
  ylabel('S : g/cm/ms')
  legend('sx','sy','sz');
  xlim(tt([1 end]))
  ylim([-18 18])
end

%%
end