function [b0map,im] = measureb0(PFILES,N,nz)
% function b0map = measureb0(PFILES,[OPTE])
% Example: b0map = measureb0(['P60416.7';'P61440.7']);
%
% Gets TE from usercv Pfile header (see
% projects/ssfpbanding/epic/ssfpbanding4.e)
%
% $Id: measureb0.m,v 1.1 2015/05/22 21:39:24 jfnielse Exp $

[im1,d] = loadim(PFILES(1,:),N,nz);
[im2,d] = loadim(PFILES(2,:),N,nz);

% slice = 5;
% im1 = im1(:,:,slice);
% im2 = im2(:,:,slice);

% hsize = 3; sigma= 1.5;
% for ii = 1:size(im1,3)
% im1(:,:,ii) = smoothim(im1(:,:,ii),hsize,sigma);
% im2(:,:,ii) = smoothim(im2(:,:,ii),hsize,sigma);
% end

% OPTE = [usercv1(2) usercv2(2)] * 1e-6           % sec
OPTE = [3, 5.3] * 1e-3 %sec

b0map = angle(im2./im1)/(OPTE(2)-OPTE(1))/(2*pi);    % Hz

im = abs(im1);
im = im./max(im(:));

return;

% EOF
