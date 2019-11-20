function [Jx, Jy, Jz, J] = formJ_ls(A, rf, mask)
% Inspired from codes inside Hao's kb_Jacobian_inputA()
% For LS objective ||sum_{coils}{SAb} - d||_2^2, Calc J == S*X*A*b W/ mask
% Basis matrix H not included as it's a degree of freedom in pulse design.
% Ref. Sun et al. DOI: 10.1109/TMI.2015.2478880

xyz_n = mask2mlocn(mask);
xyz_mn = xyz_n(~~xyz_n(:,1), 2:end); % _m: masked; _n: normalized

if isnumeric(A)
  rf = spCompatible(rf);
  J_tmp = bsxfun(@times, permute(xyz_mn, [1,3,2]), bsxfun(@times, A, rf(:).'));
  [Jx, Jy, Jz] = deal(J_tmp(:,:,1), J_tmp(:,:,2), J_tmp(:,:,3));
  J = J_tmp(:,:);
else
  arg.A = A;
  arg.rf = rf(:);
  arg.locm = xyz_mn(:,1); % for Jx
  
  Jx = fatrix2('idim',numel(rf), 'omask',mask, 'arg',arg, 'does_many',true ...
    , 'forw',@forw, 'back',@back);
  Jy = Jx; Jy.arg.locm = xyz_mn(:,2);
  Jz = Jx; Jz.arg.locm = xyz_mn(:,3);
  J = [Jx, Jy, Jz];
end

end

function y = forw(arg, x)
[locm, A, rf] = deal(arg.locm, arg.A, arg.rf);

rfx = full(bsxfun(@times, spCompatible(rf), x));
y = bsxfun(@times, locm, A*rfx);
y = embed(y, A.omask);
end

function x = back(arg, y)
[locm, A, rf] = deal(arg.locm, arg.A, arg.rf);

if size(y,1) ~= nnz(A.omask)
  if size(y,1) ~= numel(A.omask), y = reshape(y, numel(A.omask), []); end
  y = y(A.omask,:);
end

locmy = bsxfun(@times, locm, y);
x = bsxfun(@times, A'*locmy, conj(rf));

end
