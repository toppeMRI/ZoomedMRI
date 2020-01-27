function gcrush = makecrusher(ncycles,opslthick)
% 
% INPUTS:
%   ncycles     -- number of cycles of phase across slab
%   opslthick   -- cm
% 
% $Id: makecrusher.m,v 1.1 2015/05/22 21:39:23 jfnielse Exp $

gamma = 4257.5;                              % Hz/Gauss
mxg = 3.9;                                   % G/cm
mxs = 14e3;                                  % G/cm/sec
area = ncycles/(gamma*opslthick);            % G/cm*sec
gcrush = trapwave(area,4e-6,mxg,mxs);
gcrush = makeevenlength(gcrush(:));

% EOF
