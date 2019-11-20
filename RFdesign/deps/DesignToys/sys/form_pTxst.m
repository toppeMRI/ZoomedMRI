function A = form_pTxst(F, sMap, mask)
% function A = form_pTxst(F, sMap, mask)
% small tip pTx system matrix
%INPUTS:
% - F (np, nk) form_FTst matrix/fatrix2 object
%OPTIONAL:
% - sMap (nx, ny, nz, nc), coil transmit sensitivity map
% - mask (nx, ny, nz), boolean, target mask

if nargin==0, test(); return; end

if ~nonEmptyVarChk('sMap') || all(sMap(:)==1), A = F; return; end
if isscalar(sMap), warning('scalar sMap?'); A = sMap*F; return; end

nc = size(sMap, 4);
[~, nk] = size(F);
sMap = reshape(sMap, [], nc);
sMap_m = sMap(mask, :);

if isnumeric(F)
  A = bsxfun(@times, permute(sMap_m, [1,3,2]), F);
  A = A(:,:);
else
  arg.F = F;
  arg.sMap_m = sMap_m; % _m: masked
  
  % modify mex nufft before switching to 'does_many'==true
  A = fatrix2('arg',arg, 'idim',nk*nc, 'omask',mask,'does_many',false ...
    ,'forw',@forw, 'back',@back);
end
end

function y = forw(arg, x)
[F, sMap_m] = deal(arg.F, arg.sMap_m);
[ni_m, nc] = size(sMap_m);

if size(x,2) ~= nc, x = reshape(x, [], nc); end

ytmp = zeros(ni_m, nc);
if nc == 1, ytmp(:,1) = F*x(:,1);
else, parfor ic = 1:nc, ytmp(:,ic) = F*x(:,ic); end
end
y = sum(bsxfun(@times, ytmp, sMap_m), 2);

end

function x = back(arg, y)
[F, sMap_mc] = deal(arg.F, conj(arg.sMap_m));
[~, nc] = size(sMap_mc);

ytmp = bsxfun(@times, sMap_mc, y);
x = zeros(size(F,2), nc);
if nc == 1, x(:,1) = F'*ytmp(:,1);
else, parfor ic = 1:nc, x(:,ic) = F'*ytmp(:,ic); end
end
x = x(:);
end

function test()
nx = 128;
imSize = [nx,nx];
P0 = phantom(nx);
isFatrix = true; % one can also try setting this to false

[xx, yy] = deal((1:nx)/nx);
[xx, yy] = deal(xx(:), yy(:)');
sMap = exp(-sqrt(xx.^2+yy.^2)) .* exp(1i*pi*32*(xx-yy)/nx);

imMask = true(imSize);

kMask = true(imSize);
kTraj = mask2mloc(kMask); % cycle/FOV, indexed by matlab dflt order
kTraj = kTraj(:,2:end);

bP = -1i*mfft2(P0)/prod(imSize); % pulse to generate the P

% no need to concern about the Gmri warning
F_st = form_FTst(kTraj, imMask, [], 'isFatrix',isFatrix);
FS_st = form_pTxst(F_st, sMap, imMask);

P1 = reshape(FS_st*bP(:), imSize);
P2 = sMap.*reshape(F_st*bP(:), imSize);

figure, % error of order 1e-7 is due to single precision, etc.
subplot(121), imagesc(real(P1) - real(P2)); colorbar;
subplot(122), imagesc(imag(P1) - imag(P2)); colorbar;

disp('form_pTxst.test() done, spectral pulse not tested');
end
