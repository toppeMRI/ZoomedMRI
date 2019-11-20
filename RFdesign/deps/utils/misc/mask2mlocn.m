function mlocn = mask2mlocn(mask)

dims = size(mask);
if iscolumn(mask), dims = dims(1); end
mlocn = mask2mloc(mask);
mlocn(:, 2:end) = bsxfun(@times, 1./dims, mlocn(:,2:end));

end
