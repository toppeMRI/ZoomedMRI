function [e, As] = eval_func(A, b, d, sens, mask)


ncoils = size(sens,4); 
Nt = size(A,2); 
% combine sens and A into As
Atmp = [];
for ic = 1:ncoils
   sensi = sens(:,:,:,ic);
   Atmp = [Atmp; spdiag(sensi(mask))*A];
end
As = reshape(Atmp, [sum(mask(:)),ncoils*Nt]);

e = d(mask) - As*b;
