function maxMask = mask_FOV(d)
% create the mask to analysis the maximum error in the non-excited region
maxMask = zeros(size(d));
maxMask(d~=0) = 1;
ind = find(maxMask(32,32,:)~=0);
if ind(1)>1
    maxMask(:,:,ind(1)-1) = maxMask(:,:,ind(1));
end
if ind(end)<size(maxMask,3)
    maxMask(:,:,ind(end) +1) = maxMask(:,:,ind(end));
end
for imask = 1:size(maxMask,3)
    SE = ones(4,4);
    maxMask(:,:,imask) = imdilate(maxMask(:,:,imask),SE);
end
maxMask = logical(maxMask);