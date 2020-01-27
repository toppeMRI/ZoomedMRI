function g = mybridged(area_target,gbeg,gend)
% function g = mybridged(area_target,gbeg,gend)
%
% Make bridged gradient waveform.
%
% INPUTS:
%   area   - G/cm*sec
%   gbeg   - beginning value (G/cm)
%   gend   - ending value (G/cm)

mxg = 3.9;           % G/cm
mxs = 14900;       % G/cm/sec
dt  = 4e-6;        % sample duration (sec)

s = mxs * dt      % max change in g per sample (G/cm)

gpre = gbeg:s:gend;
g = gpre;
area = sum(g)*dt;
while area_target-area > 0
	gpre = [gpre min(gpre(end)+s,mxg)];
	g = [gpre gpre(end):-s:gend];
	area = sum(g)*dt;
end

return;

