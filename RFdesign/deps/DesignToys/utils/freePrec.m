function Mo = freePrec(Mi, dur, T1, T2, b0Map)
% ROTATING FRAME: decay during free precession
% The function assumes spin magnitude of 1
%INPUTS:
% - Mi (nSpins, xyz): input spins' magnetizations, w/ equilibrium norm 1s
% - T1 & T2 Sec: globally or spin-wisely.
%     (1,), global
%     (nSpins, 1), spin-wise
% - dur (1,) Sec: duration of free precesion
% - b0Map (1,) or (nSpins, 1) Hz.
%OUTPUTS:
% - Mo (nSpins, xyz)

if isequal(dur, 0), Mo = Mi; return; end
[E1, E2] = deal(dur./T1, dur./T2);

tcoeff  = exp(-[E2, E2, E1]);
Mo = bsxfun(@times, Mi, tcoeff);
Mo(:,3) = Mo(:,3) + 1-tcoeff(:,3);

if ~nonEmptyVarChk('b0Map', 'var'), return; end

Mxy = (Mo(:,1) + 1i*Mo(:,2)) .* exp(-1i*2*pi*b0Map(:)*dur); % `-1i`: MxB
[Mo(:,1), Mo(:,2)] = deal(real(Mxy), imag(Mxy));

end
