function s = g2s(g, dt)
% Slew Rate from Gradient
% INPUT
% - g (Nt, Nd) Gauss/cm
% - dt,  (1)  , optional, Sec
% OUTPUT
% - s (Nt, Nd), G/cm/Sec

dt0 = envMR('get','dt');
if ~nonEmptyVarChk('dt'), dt = dt0;  end % Sec

% g = [g; g(end,:)];
% s = diff(g, 1,1)/dt; % Gauss/cm/Sec

s = [g(1,:); diff(g, 1,1)]/dt; % Gauss/cm/Sec

end
