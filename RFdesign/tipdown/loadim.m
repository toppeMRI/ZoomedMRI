function [im,d,usercv] = loadim(pfile,N,nz,echo)
% load 2d or 3d image acquired with stfr2.e
% $Id: loadim.m,v 1.1 2015/05/22 21:39:22 jfnielse Exp $

if ~exist('echo','var')
	echo = 1;
end

coil  = 1;

% we don't use loaddab slice 0
slices = 1+(1:nz); % these are the slices we care about
rhnslices = nz+1;  % number of loaddab 'slices'

[d] = rawload_jfn(pfile,N,N,1,rhnslices,1,1,1,slices,echo);

if nz>1
	im = fftshift(ifftn(fftshift(d)));
else
	im = fftshift(ifft2(fftshift(d)));
end

return;

% EOF
