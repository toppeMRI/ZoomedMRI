% check non-convexity
% plot the cost function along one dimension of k-space to see how
% non-convex it is. 
nDD = 30; 
DD = linspace(-0.005,0.005,nDD);
ind = 100; 
[e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);    
dk = 2*(2*pi)*real(1i*J'*W*e); 
for i = 1:nDD
   i
dd = DD(i); 
k1 = k;  
k1(ind) = k1(ind) + dd; 
% k1 = k1 + dd*reshape(dk,[size(k1,1),3]); 
[e1] = eval_func(k1, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);

W = diag_sp(weightIm(mask));
cost1 = e1'*W*e1; % + R(b) term

costAll(i) = cost1; 

end
figure,plot(DD(1:nDD),costAll(1:nDD)); 
