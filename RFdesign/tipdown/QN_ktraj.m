%% optimize kt-points using quasi newton methods coded by Joel W. LeBlanc
function [k,cost] = QN_ktraj(k, d, sens, mask, weightIm)
b = 0; 
[dkx, dky, dkz, db, r, A] = kb_grad(k, b, d, sens, mask,weightIm); 
b = pinv(A)*d(mask); 
save QNtmp sens mask weightIm d; 
x0 = [k(:); real(b); imag(b)]; 

maxEval = 2000; 
numVar = 5; 
opt = LBFGSOptions;
opt.maxEval = maxEval;
opt.maxIter = 30;
opt.verbose = true;
opt.relGradTol = 0;
opt.gradTol = 0;
[x,f,g,exitFlag,output] = LBFGS(@(x) objectFun(x,sens,mask,weightIm,d),x0,opt); 
for it = 1:opt.maxIter
   cost(it) = output.iter(it).f; 
end
cost = [output.f0 cost]; 

% opt = optimset('GradObj','on');
% opt.Display = 'iter';
% [X,FVAL,EXITFLAG,OUTPUT] = fminunc(@objectFun,x0, opt);

end

function [f,g,e] = objectFun(x,sens,mask,weightIm,d)
% load QNtmp; 
kx = x(1:end/5); 
ky = x(end/5+1:end/5*2); 
kz = x(end/5*2+1:end/5*3); 
br = x(end/5*3+1:end/5*4);
bi = x(end/5*4+1:end); 
b = br + 1i*bi; 
[dkx, dky, dkz, db, e, A] = kb_grad([kx,ky,kz],b,d,sens,mask,weightIm); 
dbr = real(db); 
dbi = imag(db); 
g = [dkx; dky; dkz; dbr; dbi];
f = norm(e); 


end 