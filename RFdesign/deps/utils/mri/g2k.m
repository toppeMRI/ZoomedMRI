function kpcm = g2k(TxRx, g, dt, gam)
% function kpcm = g2k(TxRx, g, dt, gam)
% kTraj from Gradient, modified from Hao's g2k_hao().
% INPUT
%  TxRx, str,
%    Tx, k is assumed to end at the 0
%    Rx, k is assumed to start at the 0
%  g,   (Nt, Nd), G/cm
%  fov, (Nd,), cm
%  dt,  (1)  , optional, Sec
%  gam, (1)  , optional, Hz/Gauss
% OUTPUT
%  kpcm,(Nt, Nd), cycle/cm

[dt0, gam0] = envMR('get', 'dt','gam');
if ~nonEmptyVarChk('dt'),  dt  = dt0;  end % Sec
if ~nonEmptyVarChk('gam'), gam = gam0; end % Hz/Gauss

kpcm = cumsum(g) * dt * gam;
switch lower(TxRx)
  case 'tx', kpcm = reshape(bsxfun(@minus, kpcm(:,:), kpcm(end,:)), size(kpcm));
  case 'rx'  % nothing
  otherwise
    error('Transmit (tx) or Receive (rx)???');
end

end
