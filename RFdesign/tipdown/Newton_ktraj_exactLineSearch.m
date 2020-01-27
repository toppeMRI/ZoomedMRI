%% PSD (alternative update)
% slow
function [k,cost] = PSD_ktraj(k, d, sens, mask, weightIm)
b = 0;
alp = 0.0001;
alp = 0.00001;
[dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
b = pinv(A)*d(mask);

printf('iter | e | alp | norm(g)')
% for iter = 1:50
%    [dkx, dky, dkz, db, e] = kb_grad(k, b, d, sens, mask,weightIm);
%    [Hx, Hy, Hz] = getHessian(k,b,d,A,e, sens,mask,weightIm);
%    b = pinv(A)*d(mask);
%    [alp] = fminbnd(@(alp) objectFun(alp, k, b, sens,mask,weightIm,d,[dkx,dky,dkz],db,Hx, Hy, Hz),0,1,optimset('TolX',1e-4,'Display','off'));
%       norm(Hx(:))
%    alp = 1; 
%    k(:,1) = k(:,1) - alp * (Hx\dkx);
%    k(:,2) = k(:,2) - alp * (Hy\dky);
%    k(:,3) = k(:,3) - alp * (Hz\dkz);
% %    b = b -alp * db;
%    cost(iter) = norm(e);
%    printf('%d | %.5g | %.5g | %.3g \n', iter, norm(e), alp, norm([dkx; dky; dkz; db]))
% end
for iter = 1:50
   % get the Hessian matrix
   [Hx, Hy, Hz] = getHessian(k,b,d,A,e, sens,mask,weightIm); 
   b = pinv(A)*d(mask);
   [alp] = fminbnd(@(alp) objectFun(alp, k, b, sens,mask,weightIm,d,[dkx,dky,dkz],db,Hx, Hy, Hz),0,1,optimset('TolX',1e-4,'Display','off'));

   % backtracking step size
   norm(Hx(:))
   k(:,1) = k(:,1) - alp * (Hx\dkx);
   k(:,2) = k(:,2) - alp * (Hy\dky); 
   k(:,3) = k(:,3) - alp * (Hz\dkz); 
   % get the gradient 
   [dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm);
   grad_norm = norm([dkx; dky; dkz; db]); 
   alp = alp/2;
   cost(iter) = norm(e); 
   printf('%d | %.5g | %.5g | %.3g \n', iter, cost(iter), alp, grad_norm) 
end 

end

function [f] = objectFun(alp, k, b, sens,mask,weightIm,d, dk, db, Hx, Hy, Hz)
dkx = dk(:,1);
dky = dk(:,2);
dkz = dk(:,3);
k(:,1) = k(:,1) - alp * (Hx\dkx);
k(:,2) = k(:,2) - alp * (Hy\dky);
k(:,3) = k(:,3) - alp * (Hz\dkz);
% b = b-alp*db;
[dkx, dky, dkz, db, e, A] = kb_grad(k,b,d,sens,mask,weightIm);
% g = [dkx; dky; dkz; dbr; dbi];
f = norm(e);
end




function [Hx, Hy, Hz] = getHessian(k,b,d,A,e,sens,mask,weightIm)
%% construct matrices
nx = size(d,1);
ny = size(d,2);
nz = size(d,3);
x = [ceil(-nx/2):1:floor((nx-1)/2)];
y = [ceil(-ny/2):1:floor((ny-1)/2)];
z = [ceil(-nz/2):1:floor((nz-1)/2)];
[xx,yy,zz]=ndgrid(x,y,z);

kx = k(:,1); 
ky = k(:,2); 
kz = k(:,3); 

Db = diag_sp(b);
X = diag_sp(xx(mask));
Y = diag_sp(yy(mask)); 
Z = diag_sp(zz(mask)); 
Sr_tmp = sens(:,:,:,1); 
Sr = diag_sp(Sr_tmp(mask)');
W = diag_sp(weightIm(mask));

%% construct Hessians
Hmiddle = A'*(X*Sr'*W*Sr*X)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
H2 = 2*real(conj(b).*diag(Hmiddle).*b)-2*real(full(e'*W*Sr*X*X*A*Db)');
Hx = triu(H1,1)+tril(H1,-1)+diag(H2); 

Hmiddle = A'*(Y*Sr'*W*Sr*Y)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
H2 = 2*real(conj(b).*diag(Hmiddle).*b)-2*real(full(e'*W*Sr*Y*Y*A*Db)');
Hy = triu(H1,1)+tril(H1,-1)+diag(H2); 

Hmiddle = A'*(Z*Sr'*W*Sr*Z)*A; 
H1 = 2*real(full(Db'*Hmiddle*Db));
H2 = 2*real(conj(b).*diag(Hmiddle).*b)-2*real(full(e'*W*Sr*Z*Z*A*Db)');
Hz = triu(H1,1)+tril(H1,-1)+diag(H2); 

end
