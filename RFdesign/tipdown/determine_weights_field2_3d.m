function [x_est,projection,nrmse] = determine_weights_field2_3d(sx,sy,sz,kx,ky,kz,d,mask,sens_proj,fieldmap,pulse_params,lambda, miniMax,ishard)

ncoils = size(sens_proj,4);

dx = size(d,1);
dy = size(d,2);
dz = size(d,3); 

freqs.kx = kx;
freqs.ky = ky;
freqs.kz = kz;

gam = 4257*2*pi; % rad/s
dt = 4e-6; %s
%A =ExcitationMatrixField2D(d,sens_proj,mask,freqs,fieldmap,sx,sy,pulse_params);
A = [];
if ishard == 1
    durations = setup_duration_hard(freqs,sx,sy,sz,pulse_params);
else 
    durations = setup_duration(freqs,sx,sy,pulse_params);
end 
x = [0:1:dx-1]-floor(dx/2);
y = [0:1:dy-1]-floor(dy/2);
z = [0:1:dz-1]-floor(dz/2);
[xx,yy,zz] = ndgrid(x,y,z);

for i = 1:length(kx)
    w = exp(-1i*2*pi*fieldmap*durations(i));
    freq = exp(1i*2*pi*xx*kx(i)/dx+1i*2*pi*yy*ky(i)/dy+1i*2*pi*zz*kz(i)/dz);
    freq_w = freq.*w;
    for j = 1:ncoils
        a_ij = freq_w.*squeeze(sens_proj(:,:,:,j));
        A = [A a_ij(mask)];        
    end
end


d_red = d(mask);
if miniMax == 0
    x_init_inv = inv(A'*A+lambda*eye(ncoils*length(kx)))*(A'*d_red);
    nrmse_inv = norm(d_red-A*x_init_inv)/norm(d_red);
else
    [x_init_inv, nrmse_inv] = minInfiniteNorm(A,d_red);
end
%{
%n_iter = length(kx)*ncoils*2;
n_iter = 300;
x_init_iter = zeros(length(kx)*ncoils,1);
[x_init_iter,info] = qpwls_pcg1(x_init_iter,A,1,d_red,0,'niter',n_iter);   

nrmse_iter = norm(d_red-A*x_init_iter)/norm(d_red);

projection = embed(A*x_init_iter,mask);
error_vec = (d-projection);
nrmse = norm(error_vec(:))/norm(d(:));
%}
projection = embed(A*x_init_inv,mask);
nrmse = nrmse_inv;
x_est = [];
for i = 1:ncoils
    weight_i = x_init_inv(i:ncoils:end);
    x_est = [x_est weight_i(:)];
end




