function gcrush = makeCrusher(kpcm, gMax, sMax, dt, gam) 
%INPUTS:
% - kpcm (1,), cycle/cm
%OUTPUT:
% - gcrush(ng, 1), G/cm

[gMax0, sMax0, dt0, gam0] = envMR('get', 'gMax','sMax','dt','gam');
if ~nonEmptyVarChk('gMax'), gMax = gMax0; end % Gauss/cm
if ~nonEmptyVarChk('sMax'), sMax = sMax0; end % Gauss/cm/Sec
if ~nonEmptyVarChk('dt'),   dt   = dt0;   end % Sec
if ~nonEmptyVarChk('gam'),  gam  = gam0;  end % Hz/Gauss

area = kpcm/gam; % G/cm*sec
gcrush = makeTrapezoid(area,gMax,sMax,dt);

end


