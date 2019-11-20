function g = makeTrapezoid(area, gMax, sMax, dt)
%INPUTS
% - area (1,), G/cm * Sec
%OUTPUTS
% - g (ng, 1), G/cm, trapezoid (or triangle), with given area

[gMax0, sMax0, dt0] = envMR('get', 'gMax', 'sMax', 'dt');
if ~nonEmptyVarChk('gMax'), gMax = gMax0; end % Gauss/cm
if ~nonEmptyVarChk('sMax'), sMax = sMax0; end % Gauss/cm/Sec
if ~nonEmptyVarChk('dt'),   dt   = dt0;   end % Sec

a_th = gMax.^2/sMax; % _th: threshold

area_abs = abs(area);
if area_abs <= a_th
  g_p_tmp = sqrt(area_abs*sMax); % _p: peak
  n_h = ceil(g_p_tmp/sMax/dt); % _h: half
  g_p = area / n_h / dt;
  amp = [1/n_h : 1/n_h : 1, (1-1/n_h) : -1/n_h : 1/n_h];
  g = g_p * amp;
else
  area_rect = area_abs - a_th;
  n_hr = ceil(gMax/sMax/dt);
  g_pr = a_th/2 / n_hr / dt; % _r: ramp
  amp_r = (1/n_hr : 1/n_hr : 1);
  
  n_h = ceil(area_rect/gMax/dt);
  g_p = area_rect / n_h / dt;
  
  g = [g_pr * amp_r, g_p*ones(1,n_h), g_pr*amp_r(end:-1:1)];
  
end

end
