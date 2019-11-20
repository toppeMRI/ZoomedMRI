function g = k2g(TxRx, kpcm, dt, gam)
% Gradient from kTraj, modified from Hao's k2g_hao().
% INPUT
% - TxRx, str,
%     Tx, k is assumed to end at the 0, output g(end,:) == 0
%     Rx, k is assumed to start at the 0, output g(1,:) == 0
% - kpcm (Nt, Nd) cycle/cm
%OPTIONAL
% - dt,  (1,), Sec
% - gam, (1,), Hz/G
% OUTPUT
% - g (Nt, Nd), G/cm

[dt0, gam0] = envMR('get', 'dt','gam');
if ~nonEmptyVarChk('dt'),  dt  = dt0;  end % Sec
if ~nonEmptyVarChk('gam'), gam = gam0; end % Hz/Gauss

g = reshape([kpcm(1,:); diff(kpcm(:,:),1,1)]/gam/dt, size(kpcm));

switch lower(TxRx)
  case 'tx', assert(all(kpcm(end,:) == 0), 'Transmit k-space must end at 0');
  case 'rx' % nothing
  otherwise
    error('Transmit (tx) or Receive (rx)???');
end

end
