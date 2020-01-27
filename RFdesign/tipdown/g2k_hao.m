function k = g2k_hao(g)
%function k = g2k_hao(g)
% g: Nt x 3
dt = 4e-6; %sec
gamma = 4257; 
k = -flipud(cumsum(flipud(g),1)*dt*gamma); 
end