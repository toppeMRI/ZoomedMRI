function rf = rf_qpwls(phm, kpcm, d, varargin)
% solve ||Ab-d||_2^2 + beta ||b||_2^2
%INPUTS:
% - phm struct:
%   .sMap (nx, ny, nz, nc) or [], transmission sensitivity map
%   .w (nx, ny, nz) or [], spatial weightfor qpwls
%   .b0Map (nx, ny, nz) or [], Hz, b0 field map,
%   .fov (xyz,), cm
%   .ofst (xyz,)
% - kpcm (nk, xyz), cycle/cm
% - d (nx, ny, nz)
%OPTIONAL:
% - beta (1,)
% - nter (1,)
% - rf0 (nk,)
% - dt (1,)
% - gam (1,)
% - isFatrix [t/F]
%OUTPUTS:
% - rf
[sMap, w, b0Map, fov, ofst] = getattrs(phm, {'sMap','w','b0Map','fov','ofst'});

arg = envMR('get_s', 'dt', 'gam');
arg.isFatrix = false;
[arg.rf0, arg.beta, arg.niter] = deal([], 0, 20);

arg = attrParser(arg, varargin);

[rf0, beta, niter] = getattrs(arg, {'rf0','beta','niter'});
FT_c = getattrs(arg, {'dt','gam','isFatrix'}, true);

m = ~~w;
k = bsxfun(@times, kpcm, fov); % cycle/cm -> cycle/FOV

F = form_FTst(k, m, b0Map, 'ofst',ofst,FT_c{:});
A = form_pTxst(F, sMap);

if isempty(m), m = ':'; end % trick, x(':') == x(:)

W = Gdiag(w(m));
rf = qpwls_pcg1(rf0, A, W, d(m), sqrt(beta), 'niter',niter);

% res = embed(A*rf, m); figure, im(res); % test
end

