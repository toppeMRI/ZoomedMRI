function F = form_FTst(k_t, mask, b0Map_b0rng, varargin)
% FTst = form_FTst(k_t, mask, b0Map_b0rng, dt, gam, isFatrix)
% Basically a Gmri adjoint wrap
% Purely based on small tip (ST) assumption, that all m0 are aligned w/ z-axis
% at the very beginning.
% It's unlikely that input k-space is Cartesian, thus not optimized for it
%INPUTS
% - k_t, k-traj or pulse length
%    k ~ kspace (nk, nd): cycle/fov, discrete description of small-tip k-space
%     trajectory. If plotted against time along dim1, these k-space points are
%     seperated by dt
%    t ~ trf (1,): Sec, for spectral pulse
% - mask (nx, ny, nz) or []: boolean; input [] if for spectral pulse
% - b0Map_b0rng, b0Map or off-resonance range (relaxation is not considered).
%    b0Map (nx, ny, nz): Hz, off-reconance map
%    b0rng (nfreq,): Hz, when k_t ~ t, for spectral pulse, e.g. b0rng = -30:29;
%OPTIONAL
% - ofst (xyz,), a.u. [-0.5, 0.5], dflt: 0;
% - isFatrix [t/F]
% - dt (1,),  Sec, the delta that matches the dt of kspace
% - gam (1,), Hz/G, gyro frequency
%OUTPUTS
% - F (np, nk), matrix/fatrix2 object, system matrix of small tip pulse design

if nargin == 0, test(); return; end

arg = envMR('get_s', 'dt','gam');
arg.ofst = 0;
arg.isFatrix = false;
arg = attrParser(arg, varargin);

[dt, gam, ofst, isFatrix] = getattrs(arg, {'dt','gam','ofst','isFatrix'});

gambar = 2*pi * gam;
if isscalar(k_t) % prep arg for pure spectral pulse system matrix
  trf = k_t;
  mask = true(numel(b0Map_b0rng), 1); % b0rng
  ti = -trf:dt:0;
  Db0 = max(b0Map_b0rng) - min(b0Map_b0rng);
  warning('Feature not tested, constant coeff may be incorrect');
  expb0shift = exp(1i*2*pi*(min(b0Map_b0rng) - ctrInd(Db0)));
  
  k = Db0 * ti.';
  [b0Map, ti] = deal([]);
else
  trf = dt*(size(k_t,1)-1);
  k = k_t;
  if ~nonEmptyVarChk('b0Map_b0rng') || ~any(b0Map_b0rng(:))
    [b0Map,ti] = deal([]);
  else, [b0Map, ti] = deal(b0Map_b0rng,  (-trf:dt:0)); % think about M x B
  end
  expb0shift = 1;
end

% Effectively the DFT matrix formed becomes its own conj, fitting ST model.
% k = -k;

if ~isFatrix
  xyz_n = mask2mlocn(mask); % unit: cycle/FOV; _n: normalized.
  xyz_mn = bsxfun(@plus, xyz_n(~~xyz_n(:,1), 2:end), ofst(:)'); % _m: masked.
  if isempty(b0Map), freq0 = 0;
  else, freq0 = b0Map(mask)*ti(:)';
  end
  freq = xyz_mn * k';
  F = exp(1i*2*pi * (freq0 + freq));
else
  error("isFatrix==true bug not fixed")
  L = 6;
  Nd = size(mask);
  if iscolumn(mask), Nd = Nd(1); end
  
  Jd = 6*ones(size(Nd));
  n_shift = Nd/2; % not sure if n_shift fails for odd Nd
  nufftArgs = {Nd, Jd, 2*Nd, n_shift, 'table', 2^10, 'minmax:kb'};
  
  F = Gmri(k, mask, 'ti',ti,'zmap',b0Map, 'L',L,'nufft',nufftArgs)';
end

F = (1i*gambar*dt*expb0shift) * F;
end

function test()
nx = 128;
imSize = [nx,nx];
P0 = phantom(nx);
imMask = true(imSize);
isFatrix = true; % one can also try setting this to false

[dt, gam] = envMR('get', 'dt','gam');
gambar = 2*pi*gam;
expb0shift = 1;

kMask = true(imSize);
kTraj = mask2mloc(kMask); % cycle/FOV, indexed by matlab dflt order
kTraj = kTraj(:,2:end);

bP = -1i*mfft2(P0)/prod(imSize); % pulse to generate the P
P2 = 1i*gambar*dt*expb0shift*mifft2(bP)*prod(imSize);

% no need to concern about the Gmri warning
F_st = form_FTst(kTraj, imMask, [],'dt',dt, 'gam',gam, 'isFatrix',isFatrix);

P1 = reshape(F_st*bP(:), imSize);

figure, % error of order 1e-7 is due to single precision of Gmri, etc.
subplot(121), imagesc(real(P1) - real(P2)); colorbar;
subplot(122), imagesc(imag(P1) - imag(P2)); colorbar;

disp('form_FTst.test() done, spectral pulse not tested');
end
