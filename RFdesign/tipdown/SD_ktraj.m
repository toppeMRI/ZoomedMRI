%% steepest gradient descent
function [k,cost] = GD_steepest_ktraj(k, d, sens, mask, weightIm)
b = 0; 
alp = 0.0001; 
alp = 0.00001; 
[dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm); 
b = pinv(A)*d(mask); 

   printf('iter | e | alp | norm(g)')
for iter = 1:100
   [dkx, dky, dkz, db, e] = kb_grad(k, b, d, sens, mask,weightIm); 
   [alp fval, exitflag] = fminbnd(@(alp) objectFun(alp, k, b, sens,mask,weightIm,d,[dkx,dky,dkz],db),0,2.5e-4,optimset('TolX',1e-6,'Display','off')); 
   k(:,1) = k(:,1) - alp * dkx;
   k(:,2) = k(:,2) - alp * dky; 
   k(:,3) = k(:,3) - alp * dkz; 
   b = b -alp * db; 
   cost(iter) = norm(e); 
   printf('%d | %.5g | %.5g | %.3g \n', iter, norm(e), alp, norm([dkx; dky; dkz; db])) 
end 
 
end 

function [f] = objectFun(alp, k, b, sens,mask,weightIm,d, dk, db)
dkx = dk(:,1); 
dky = dk(:,2); 
dkz = dk(:,3); 
k(:,1) = k(:,1) - alp * dkx;
k(:,2) = k(:,2) - alp * dky;
k(:,3) = k(:,3) - alp * dkz;
b = b-alp*db; 
[dkx, dky, dkz, db, e, A] = kb_grad(k,b,d,sens,mask,weightIm);
% g = [dkx; dky; dkz; dbr; dbi];
f = norm(e);


end 

function testLineSearch(alp, k, b, sens,mask,weightIm,d, dk, db)
alpArray = 0:0.00001:0.001; 
for ia = 1:length(alpArray)
   alp = alpArray(ia);
   f(ia) = objectFun(alp, k, b, sens,mask,weightIm,d, [dkx,dky,dkz], db); 
end
figure,plot(alpArray, f); 
end 