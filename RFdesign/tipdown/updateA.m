function [A] = updateA(A,dk,mask,xx,yy,zz,nx,ny,nz)
% Update the system matrix due to k update.
% INPUTS:
%     A = system matrix (with or without fieldmap)
%     dk = k-space trajectory change
%     see formA.m

A = A.*exp(1i*2*pi*(xx(mask(:))*(dk(:,1)'/nx) + yy(mask(:))*(dk(:,2)'/ny) + zz(mask(:))*(dk(:,3)'/nz)));
