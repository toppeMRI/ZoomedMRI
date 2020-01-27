function gcrush = makecrusher_hao(ncycles,opslthick,prearea)
% 
% INPUTS:
%   ncycles     -- number of cycles of phase across slab
%   opslthick   -- cm
%   prearea     -- net gradient area in other part (Hao)
% $Id: makecrusher_hao.m,v 1.1 2015/05/22 21:39:23 jfnielse Exp $

if ~exist('prearea','var')
    prearea = 0; 
end 
gamma = 4257.5;                              % Hz/Gauss
mxg = 3.9;                                   % G/cm
mxs = 14e3;                                  % G/cm/sec
area = ncycles/(gamma*opslthick) - prearea;            % G/cm*sec
gcrush = trapwave(area,4e-6,mxg,mxs);
gcrush = makeevenlength(gcrush(:));

% EOF
