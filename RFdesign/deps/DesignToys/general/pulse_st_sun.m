function [rf, g, info] = pulse_st_sun(phm, d, isBalanced, initMeth, varargin)
%INPUTS
% - phm
%   .sMap
%   .w
%   .fov
%   .b0Map
%   .ofst  (xyz,) a.u., not cm
% - d

%% input parsing
arg = envMR('get_s'); % struct of {'dt', 'gam', 'gMax', 'sMax'}
arg.isFatrix = false; % existing fatrix2 is slow for pulse design at small size

[arg.eta, arg.ncycle, arg.nIP] = deal(0.01, 100, 20);
[arg.nCG] = deal(400);

arg.methPara = dflt_initPara(initMeth);

arg = attrParser(arg, varargin);
[dt, gam, gMax, sMax] = getattrs(arg, {'dt','gam','gMax','sMax'});
isFatrix = arg.isFatrix;

methPara = arg.methPara;
[eta, ncycle, nIP] = getattrs(arg, {'eta','ncycle','nIP'});
nCG = arg.nCG;

%% start design
d = d/max(d(:)); % small tip design, need to scale the pulse outside

%% k_pcm initialization
switch lower(initMeth)
  case {'xkt', 'extendedktpoints', 'xktpts'} % extended kt pts
    kpcm = kTx.init.extendedKtPts(phm, d, methPara.npe, methPara.tol...
                                  , gMax, sMax, dt);
  case {'spins'} % spins init
    kpcm = kTx.init.spins(methPara.dur, methPara.kmax, methPara.u...
                          , methPara.v, methPara.alpha, methPara.beta...
                          , methPara.doFast, gMax, sMax, dt);
  otherwise, error('Not supported init method')
end
% ko_pcm = kTx.init.spins(dur, kMax, u, v, alpha, beta, doFast, gMax, sMax, dt)
kpcm = kTx.projKpcm(kpcm, gMax, sMax, dt, gam); % project to feasible K
info.k_ini = kpcm;

%% optimize directly on the final k-traj (using basis function)
ki_pcm = kpcm;
env_Args = {'dt',dt, 'gam',gam, 'gMax',gMax, 'sMax',sMax};
bSIP_args = {'eta',eta, 'ncycle',ncycle, 'niter_o',nIP, 'isFatrix',isFatrix};

kpcm = kTx.opt.BsplInteriorPoint(phm, d, ki_pcm, bSIP_args{:}, env_Args{:});

%% project to feasible K
kpcm = kTx.projKpcm(kpcm, gMax, sMax, dt, gam);
g = k2g('tx', kpcm, dt, gam);

%% ramp up gradient so slew rate is not violated
g_ramp = makeRamp(g(1,:), true, sMax, dt);
g = [g_ramp; g];

%% make balance when requested
if isBalanced
  k_offset = sum(g,1)*dt; % cycle/cm
  mk_offset = max(abs(k_offset));
  
  g_b_tmp = makeCrusher(mk_offset, gMax, sMax, dt, gam)'; % _b: balancer
  g_b = bsxfun(@times, g_b_tmp, k_offset/mk_offset);
  g = [g_b; g];
end

%% Generate RF
kpcm = g2k('tx', g, dt, gam);
rf = rf_qpwls(phm, kpcm, d, 'beta',eta*size(kpcm,1),'niter',nCG,'rf0',[]...
              ,'dt',dt,'gam',gam);

info.k_res = kpcm;
end

%% Default init parameters
function methPara = dflt_initPara(initMeth)
% following parameteres chosen empirically by Sun
switch lower(initMeth)
  case {'xkt', 'extendedktpoints', 'xktpts'}
    [methPara.npe, methPara.tol] = deal(70, 0.2);
  case {'spins'}
    % sec, 3.8e-3 from Hao's Joint paper; Maximum extent of k-space in rad m^-1
    [methPara.dur, methPara.kmax] = deal(3.8e-3, 50*3);
    % Polar angular velocity default 8*pi % Azimuthal angular velocity dflt 2*pi
    [methPara.u,   methPara.v]    = deal(18*pi/methPara.dur, 2*pi/methPara.dur);
    % [Speed, position] of transition between slow and fast radial phase
    [methPara.alpha, methPara.beta] = deal(20, 0.25);
    [methPara.doFast] = deal(false);
  otherwise, error('Not supported init method')
end

end

