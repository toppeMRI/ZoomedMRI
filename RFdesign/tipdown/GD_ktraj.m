%% gradient descent
function [k,cost] = GD_ktraj(k, d, sens, mask, weightIm)
b = 0; 
alp = 0.0001; 
alp = 0.00001; 
[dkx, dky, dkz, db, e, A] = kb_grad(k, b, d, sens, mask,weightIm); 
b = pinv(A)*d(mask); 

   printf('iter | e | alp | norm(g)')
for iter = 1:100
   [dkx, dky, dkz, db, e] = kb_grad(k, b, d, sens, mask,weightIm); 
   k(:,1) = k(:,1) - alp * dkx;
   k(:,2) = k(:,2) - alp * dky; 
   k(:,3) = k(:,3) - alp * dkz; 
   b = b -alp * db; 
   cost(iter) = norm(e); 
   printf('%d | %.5g | %.5g | %.3g \n', iter, cost(iter), alp, norm([dkx; dky; dkz; db])) 
end 
   