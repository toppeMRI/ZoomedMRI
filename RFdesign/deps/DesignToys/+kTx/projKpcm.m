%% project the K traj to feasible region, by 2nd-order bSpline parameterization
function k_pcm = projKpcm(k_pcm, gMax, sMax, dt, gam)
% k = projKpcm(k, Smax, Gmax)
%INPUTS:
% - k_pcm (nk, ndim), cycle/cm
% - gMax (1,), Gauss/cm: gradient constraint
% - sMax (1,), Gauss/cm/s: slew rate constraint
% - dt   (1,), Sec: dwell time
% - gam  (1,), Hz/Gauss: gyro freq
%OUTPUTS:
% - k_pcm

if nargin == 0, test(); return; end

% start from here
[gMax0, sMax0, dt0, gam0] = envMR('get', 'gMax', 'sMax', 'dt', 'gam');
if ~nonEmptyVarChk('gMax'), gMax = gMax0; end % Gauss/cm
if ~nonEmptyVarChk('sMax'), sMax = sMax0; end % Gauss/cm/Sec
if ~nonEmptyVarChk('dt'),   dt   = dt0;   end % Sec
if ~nonEmptyVarChk('gam'),  gam  = gam0;  end % Hz/Gauss
if sMax < 100, sMax = sMax*1000; end % convert legacy input unit to Gauss/cm/s

nCycle = 200; % chosen empirically by Sun Hao

[nk, ndim] = size(k_pcm);
[B, dB1, dB2] = kTx.basis.Bspline2ndOrder_regular(nk, nCycle, ndim, dt, gam);
[dB, bnd] = deal([dB1;dB2],[gMax*ones(size(dB1,1),1);sMax*ones(size(dB2,1),1)]);

coeff = B\k_pcm;

% Theoretically the quadprog should have following inputs. eye() is an approx of
% BtBk = kron(speye(3), spCompatible(B'*B));
% [coeff,~,~,~] = quadprog(BtB_k, -BtB_k*coeff(:), dB, bnd);
% quadOptions = optimset('Algorithm', 'interior-point-convex'); % dflt already
[coeff,~,~,~] = quadprog(speye(numel(coeff)), -coeff(:), dB, bnd);
k_pcm = B * reshape(coeff, size(B,2), []);

end

%%
function test()

k_pcm_raw = (rand(999, 3)-0.5) * 4; % cycle/cm
k_pcm_raw = [k_pcm_raw; zeros(1,3)];

[gMax, sMax, dt, gam] = envMR('get', 'gMax', 'sMax', 'dt', 'gam');
k_pcm_proj = kTx.projKpcm(k_pcm_raw, gMax, sMax, dt, gam);

g = k2g('tx', k_pcm_proj, dt, gam);
s = g2s(g, dt);

% figure
% subplot(211), plot(g); title(['gMax: ', num2str(gMax)]);
% subplot(212), plot(s); title(['sMax: ', num2str(sMax)]);

assert(all(max(abs(g), [], 1) <= gMax), 'gMax violated.');
assert(all(max(abs(s), [], 1) <= sMax), 'sMax violated.');

end
