function [Mo, Mhst] = blochSim(Mi, Beff, T1, T2, dt, gam, doCim, doGPU)
% bloch Simulation for multiple spin with globally or spin-wisely prescripted
% parameters and pulse
%
%INPUTS:
% - Mi (nSpins, xyz): input spins' magnetizations
% - Beff, Gauss:
%     (nSpins, xyz, nSteps), spin-wise;
%     (nSteps, xyz) global
% - T1 & T2, Sec: globally or spin-wisely defined T1 and T2
%     (1,), global
%     (nSpins, 1), spin-wise
% - dt (1,) Sec: temporal simulation step size.
% - gam, Hz/G, gyro frequency ratio:
%     (1,), global
%     (nSpins, 1), spin-wise NOT FINISHED
% - doCim (T/[F]), call c-based simulation
% - doGPU (T/[F]), use matlab GPU features
%OUTPUTS:
% If doGPU is enabled, the OUTPUTs will be as gpuArray
% - Mo   (nSpins, xyz, nSteps): output spins' magnetizations
% - Mhst (nSpins, xyz, nSteps): output spins' magnetizations history
%Notices:
% 1. Not much sanity check inside this function, user is responsible for
%    matching up the dimensions.
% 2. For forking, user may want to avoid using cross() in MatLab, which involves
%    many permute()-calls and is time-consuming
% 3. Put decays at the end of each time step may still be problematic, since
%    physically the spin decays continuously, this noise/nuance may worth study
%    for applications like fingerprinting simulations, etc.

if nargin == 0, test(); return; end

[dt0, gam0] = envMR('get', 'dt', 'gam');
if ~nonEmptyVarChk('dt'),  dt  = dt0;  end % Sec
if ~nonEmptyVarChk('gam'), gam = gam0; end % Hz/G

if size(Beff,3) == 1 && size(Beff,1) ~= size(Mi,1)
  % assume no single step sim for single spin
  Beff = permute(Beff, [3,2,1]); % -> (1, xyz, nSteps)
  disp('beff not being spin-specific, assuming global');
end
if isempty(Mi) || isequal(Mi,0)
  Mi = [zeros(size(Beff,1),2), ones(size(Beff,1),1)];
end
if ~nonEmptyVarChk('doCim'), doCim = false; end
if ~nonEmptyVarChk('doGPU'), doGPU = false; end
doHist = (nargout > 1);

if doCim
  % harder to deal w/ dim in C, so deal dim in matlab before pass in
  if nargout == 1, Mo = blochCim(Mi, double(Beff), T1, T2, dt, gam);
  else,            [Mo, Mhst] = blochCim(Mi, double(Beff), T1, T2, dt, gam);
  end
  return;
end

gambar = (2*pi)*gam;
% in unit, convert relaxations into losses/recovery per step
E1 = exp(- dt ./ T1(:));
E2 = exp(- dt ./ T2(:));

nSpins = size(Mi, 1);
nStep  = size(Beff, 3);

% negative to correct cross product (x-prod) direction;
Bmag = sqrt(sum(Beff.*Beff, 2)); % self multiplication is faster than power
niBmag = -1./Bmag;
niBmag(isinf(niBmag)) = 0;
Bn   = bsxfun(@times, Beff, niBmag); % Bn, normalized Beff;

theta = bsxfun(@times, Bmag, gambar)*dt;
Ct = cos(theta); % (nSpins, 1, nSteps) or (1, 1, nSteps)
Ct_1 = 1 - Ct;   % trick for Rotation Matrix
St = sin(theta); % (nSpins, 1, nSteps) or (1, 1, nSteps)

StBn   = bsxfun(@times, Bn, St);   % trick for Rotation Matrix
CtBn_1 = bsxfun(@times, Bn, Ct_1); % trick for Rotation Matrix

[Mx0, My0, Mz0] = deal(Mi(:,1), Mi(:,2), Mi(:,3));
if doHist, [Mox, Moy, Moz] = deal(zeros(nSpins, nStep)); end

doFlag = any(niBmag ~= 0, 1);

% Updates are calculated with rotation-then-decay, instead of the canonical
% differential equation expression.
% Discretizing differential equation may cause precision issue.

if doGPU
  fn = @gpuArray;
  [E1,  E2]           = deal(fn(E1),  fn(E2));
  [Bn,  StBn, CtBn_1] = deal(fn(Bn),  fn(StBn), fn(CtBn_1));
  [Mx0, My0,  Mz0]    = deal(fn(Mx0), fn(My0),  fn(Mz0));
  if doHist, [Mox, Moy,  Moz] = deal(fn(Mox), fn(Moy), fn(Moz)); end
end

% full lower-case variables are local in loop, o.w. not local
for istep = 1:nStep
  if doFlag(istep)
    % step-wisely extract pre-processed variables
    bn     = Bn(:,:,istep);
    stbn   = StBn(:,:,istep);
    ctbn_1 = CtBn_1(:,:,istep);

    ip = sum(bsxfun(@times, bn, [Mx0, My0, Mz0]), 2); % vector inner product
    ct = Ct(:, :, istep);

    % explicitly express cross(bn, Mo_ii_1) as a matrix vector multiplication
    mx1 =  ct       .*Mx0 -stbn(:,3).*My0 +stbn(:,2).*Mz0 +ip.*ctbn_1(:,1);
    my1 =  stbn(:,3).*Mx0 +ct       .*My0 -stbn(:,1).*Mz0 +ip.*ctbn_1(:,2);
    mz1 = -stbn(:,2).*Mx0 +stbn(:,1).*My0 +ct       .*Mz0 +ip.*ctbn_1(:,3);
  else
    mx1 = Mx0;
    my1 = My0;
    mz1 = Mz0;
  end
  % relaxation effects: "1" in Mz0 since M0=1 by assumption
  % also, update Mo_ii_1;
  Mx0 = mx1.*E2;
  My0 = my1.*E2;
  Mz0 = mz1.*E1 + 1-E1;
  if doHist, [Mox(:,istep),Moy(:,istep),Moz(:,istep)] = deal(Mx0,My0,Mz0); end
end
Mo = [Mx0, My0, Mz0];
if doHist, Mhst = permute(cat(3, Mox, Moy, Moz), [1,3,2]); end

end

%%
function test()
Mi = [1,0,0; 0,1,0; 0,0,1];
nt = 512;
t = (0:nt-1)';
Beff = 10 * [cos(t/nt * 2*pi), sin(t/nt * 2*pi), atan(t-round(nt/2))/pi];

%% test -- 1
Mo_S = blochSim(Mi, Beff, 1, 4e-2, [],[], false); % T1/T2: 1000/40 ms
tmp = Mo_S - [ 0.559535641648385,  0.663342640621335,  0.416341441715101;
               0.391994737048090,  0.210182892388552, -0.860954821972489;
              -0.677062008711222,  0.673391604920576, -0.143262993311057];


assert(norm(tmp(:)./Mo_S(:))<=1e-9);
[~, Mhst] = blochSim(Mi, Beff, 1, 4e-2, [],[], false); % T1/T2: 1000/40 ms
tmp = Mo_S - Mhst(:,:,end);
assert(norm(tmp(:)./Mo_S(:))<=1e-9);

Mo_C = blochSim(Mi, Beff, 1, 4e-2, [],[], true); % T1/T2: 1000/40 ms
[~, Mhst] = blochSim(Mi, Beff, 1, 4e-2, [],[], true); % T1/T2: 1000/40 ms
tmp = Mo_C - Mhst(:,:,end);
assert(norm(tmp(:)./Mo_C(:))<=1e-9);
tmp = Mo_S - Mo_C;
assert(norm(tmp(:)./Mo_S(:))<=1e-9);

%%
disp([mfilename, '.test() done']);
% keyboard;
end
