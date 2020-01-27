function [A] = formA_GMRI(k, sens, mask,fmap,tt, xyzRange)
% construct the nufft fatrix object. Not tested to corrected yet. 
% but seems forward operation using this is slower than direct matrix
% matrix vector multiplication

   dt = 4e-6; % seconds, RF and gradient sampling period

   Nc = size(sens,4); 
   
   dim = length(xyzRange.x); % dimension of square x-y grid
   dimz = length(xyzRange.z);
   FOV = -xyzRange.x(1)*2;
   FOVz = xyzRange.z(end) - xyzRange.z(1);
   
   % nufft params
   J = 6;K = 2*dim;Kz = 2*dimz; Jz = 6; L = 4;
   %nufft_args = {[dim dim],[J J],[K K],[dim dim]/2,'minmax:kb'};
   nufft_args = {[dim dim dimz],[J J Jz],[K K Kz],[dim dim dimz-1]/2,'minmax:kb'};
   
   gambar = 4257.6;	% gamma/2pi in Hz/g
   gam = gambar*2*pi; % gamma in radians/g
   % trick: the system matrix is just the transpose of a SENSE image recon matrix!
   %Gsml = Gmri_SENSE(k,mask,'fov',[FOV FOV],'basis',{'dirac'}, ...
   
%    save tmpCont
%    load tmpCont

   %Hao: kspace have to be rescaled here.
   sens_reshape = reshape(sens,[dim*dim*dimz Nc]);
   sens_masked = sens_reshape(mask,:);
   FOVz_s = FOVz*dimz/(dimz-1);
   %    FOVz_s = FOVz;
   A = Gmri_SENSE(k,logical(mask),'fov',[FOV FOV FOVz_s],'basis',{'dirac'}, ...
      'nufft',nufft_args,'exact',0,... % try exact = 1 here;
      'sens',conj(sens_masked)*(-1i*gam*dt), ...
      'ti',tt,'L',L,'zmap',2*pi*1i*fmap)';
