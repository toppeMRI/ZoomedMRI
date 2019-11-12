

load P51200_res.mat
P = sqrt(sum(abs(fPS(:,:,:,:,1)).^2,4));
kP2 = zeros(80,80,80);
kP = fftshift(fftn(fftshift(P)));
kP2(:,:,(end/2-20):(end/2+19)) = kP;
