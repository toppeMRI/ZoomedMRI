function gRamp = makeRamp(g0, isUp, sMax, dt)
% for the sake of hardware, ramp the gradient from 0
%INPUTS:
% - g0 (xyz,), Gauss/cm
% - isUp (T/f), ramp UP or down
% - sMax (1,), Gauss/cm/Sec
%OUTPUTS:
% - gRamp (xyz,), gradient ramp

if ~nonEmptyVarChk('isUp'), isUp = true; end

[sMax0, dt0] = envMR('get', 'sMax', 'dt');
if ~nonEmptyVarChk('sMax'), sMax = sMax0; end % Gauss/cm/Sec
if ~nonEmptyVarChk('dt'),   dt = dt0; end % Sec
sMax = min(sMax)*dt; % in case of sMax not isotropic, take minimal one

ramp_nt = ceil(max(abs(g0))./sMax)+1;
ramp_raw = linspace(0, 1-1/ramp_nt, ramp_nt-1); % no need to keep 1

gRamp = ramp_raw(:)*g0(:).';

if ~isUp, gRamp = flipud(gRamp); end % ramp down

end
