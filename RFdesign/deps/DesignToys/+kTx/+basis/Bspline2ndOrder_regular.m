function [B, dB1, dB2] = Bspline2ndOrder_regular(nk, ncycle, ndim, dt, gam)
% function [B, dB1, dB2] = Bspline2ndOrder_regular(nk, ncycle, ndim)
%INPUTS:
% - nk     (1,)
%OPTIONAL
% - ncycle (1,)
% - dt     (1,), Sec
% - gam    (1,), Hz/Gauss
%OUTPUTS:
% - B   (nk, ncycle),                      cycle/cm
% - dB1 (ndim*nXtrm1, ndim*nBasis) sparse, Gauss/cm
% - dB2 (ndim*nXtrm2, ndim*nBasis) sparse, Gauss/cm/Sec

if ~nonEmptyVarChk('ncycle'), ncycle = 100; end
if ~nonEmptyVarChk('ndim'),   ndim = 3; end

[dt0, gam0] = envMR('get', 'dt','gam');
if ~nonEmptyVarChk('dt'),  dt  = dt0;  end % Sec
if ~nonEmptyVarChk('gam'), gam = gam0; end % Hz/Gauss

nd_fn = @(dB)kron(speye(ndim), spCompatible([dB; -dB]));

tis = repmat(linspace(0,ncycle+1,nk)', [1,ncycle]);
B = spCompatible(bspline_1d_synth(eye(ncycle), tis,'ending','zero','order',2));

dB1_raw = k2g('tx', B, dt, gam); % Gauss/cm
dB2_raw = g2s(dB1_raw, dt);      % Gauss/cm/Sec

%% pick extremes and kron-prod
[~, imax] = max(dB1_raw);
[~, imin] = min(dB1_raw);
xtrmInd1 = union(imax, imin);

[~, imax] = max(dB2_raw);
[~, imin] = min(dB2_raw);
xtrmInd2 = union(imax, imin);

[dB1, dB2] = deal(nd_fn(dB1_raw(xtrmInd1,:)), nd_fn(dB2_raw(xtrmInd2,:)));
end
