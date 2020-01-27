% ADMM_ktraj.m
function [k,cost] = ADMM_ktraj(k, d, sens,mask,weightIm,beta,xfov,yfov,zfov)
N = length(k);
C = Cfd_all(N);
tic
nout = 50;
printf('iout | e')
lambda = 0; 
mu = 1; 
for iter = 1:nout
    [k,cost] = ADMM_updateK_ktraj(k, d, sens,mask,weightIm,beta,lambda,mu,C,xfov,yfov,zfov);
    lambda = max(0,lambda - mu*C*k(:)); 
    printf('%d | %.5g\n', iter, cost(end))
end
exet = toc
end


function C = Cfd_all(N)
C1 = Cfd1_zeros(N);
C2 = C1*C1;
forwc = @(arg, x) [C1*x(1:end/3);C1*x(end/3+1:end/3*2);C1*x(end/3*2+1:end); ...
    C2*x(1:end/3);C2*x(end/3+1:end/3*2);C2*x(end/3*2+1:end)];
backc = @(arg, x) [C1'*x(1:end/6)+C2'*x(end/6*3+1:end/6*4);...
    C1'*x(end/6+1:end/6*2)+C2'*x(end/6*4+1:end/6*5);...
    C1'*x(end/6*2+1:end/6*3)+C2'*x(end/6*5+1:end)];
C = fatrix2('forw', forwc, 'back', backc, 'idim', 3*N, 'odim', 6*N);
end

