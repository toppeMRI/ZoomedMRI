function [b1,gx,gy,gz,projIm] = cont_gen_rf_exactA(d,fieldmap,mask,sens,xyzRange,weightIm, k, B1_init, W)
betaV = 8; 

ncoils = size(sens,1); 
if(~isvar('B1_init'))
   B1_init = zeros(size(k,1),ncoils); 
end

dt = 4e-6; tt = [0:length(k)-1]*dt-length(k)*dt;

xfov = xyzRange.x(end) - xyzRange.x(1) + xyzRange.x(2) - xyzRange.x(1);
yfov = xyzRange.y(end) - xyzRange.y(1) + xyzRange.y(2) - xyzRange.y(1);
zfov = xyzRange.z(end) - xyzRange.z(1) + xyzRange.z(2) - xyzRange.z(1);

A = formA([k(:,1)*xfov, k(:,2)*yfov, k(:,3)*zfov], permute(sens,[2,3,4,1]), mask, fieldmap, tt, xyzRange);
% combine sens and A into As
if ncoils ==1
   As = A;
else
   Atmp = [];
   for ic = 1:ncoils
      sensi = sens(ic,:,:,:);
      Atmp = [Atmp; spdiag(sensi(mask))*A];
   end
   As = reshape(Atmp, [sum(mask(:)),ncoils*length(kx)]);
end

betavec = sqrt(betaV)*ones(size(As,2),1);
W = diag_sp(weightIm(mask)); 
b1 = qpwls_pcg1_hao(B1_init, As, W, d(mask), diag_sp(betavec),'niter',100);

g = k2g_hao(k);
gx = g(:,1);
gy = g(:,2);
gz = g(:,3);

projIm = embed(As*b1,mask); 
