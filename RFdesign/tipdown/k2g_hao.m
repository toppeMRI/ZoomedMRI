function g = k2g_hao(k)
% dimension of k is Nt by 3
dt = 4e-6; % s
gamma = 4257; % Hz/Gauss

kPadding = [k; k(end,:)]; 
g = diff(kPadding,1,1)/gamma/dt; 
end