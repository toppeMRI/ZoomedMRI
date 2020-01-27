% compare the analytical gradient to the numerical gradient
nInd = 50; 
for i = 1:nInd
ind = i; 
dd = 1e-8; 
k1 = k; 
k2 = k1; 
k2(ind) = k2(ind) + dd; 

[e1] = eval_func(k1, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
[e2] = eval_func(k2, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
W = diag_sp(weightIm(mask));
cost1 = e1'*W*e1; % + R(b) term
cost2 = e2'*W*e2; % + R(b) term 
dkNum = (cost2-cost1)/dd; 
dkNumA(i) = dkNum; 

[e, As, J] = kb_Jacobian_withB0(k, b, d, sens, mask,weightIm,beta, fieldmap,tt, xyzRange);
dk = 2*(2*pi)*real(1i*J'*W*e);

end
figure,plot(dkNumA); hold on; plot(dk(1:nInd),'r')
figure,plot(dk(1:nInd)./dkNumA(1:nInd)'); 
dk(ind)
