function [M]=parallel_blochsim_field(Minit,pulses,gx,gy,gz,sens,x,y,z,dt,fieldmap,fieldmap_start)
%set up bz with gradients
[X,Y,Z] = meshgrid(x,y,z);
X = X(:)';
Y = Y(:)';
Z = Z(:)';

gam=26751;                   %rad/sec/g
gambar = gam/2/pi;           % Hz/g
th_grad = fieldmap/gambar;
th_grad = th_grad(:)';
%bz = (gx*X+gy*Y+gz*Z+ones(size(gx))*th_grad)*10^(-4);
bz = gx*X;
bz = bz+gy*Y;
bz = bz+gz*Z;
eff = [zeros(fieldmap_start,1);ones(length(gz)-fieldmap_start,1)];
bz = bz+eff*th_grad;
bz = bz*10^(-4);


%set up bx and by with pulses and sens
bx = 0;
by = 0;
for i = 1:size(sens,3)
    pulse_i = squeeze(pulses(:,i));
    sens_i = squeeze(sens(:,:,i));
    sens_i_ext = [];
    for j = 1:length(z)
        sens_i_ext(:,:,j)=sens_i;
    end
    sens_i_ext = sens_i_ext(:);
    b_xy = pulse_i*(sens_i_ext.');
    bx = bx+real(b_xy);
    by = by+imag(b_xy);
end



T1 = 1000*ones(length(X),1);
T2 = 100*ones(length(X),1);
if (Minit == 0)
    Mi = [0;0;1]*ones(1,length(X));
else
    Mi = Minit(:)';
    Mi = [zeros(size(Mi));zeros(size(Mi));Mi];
end

[mx,my,mz] = blochsim2(Mi, bx, by, bz, T1, T2, dt*10^(3));

m_trav = mx(end,:)+1i*my(end,:);
M = reshape(m_trav, length(x),length(y),length(z));
