function demo()

% `b0Map_hi`, `dIV_hi`, and `m_hi` should have the same matrix size.
f = matfile('./b0Map.mat');  % off-resonance map
b0Map_hi = f.b0Map;

f = matfile('./iv.mat');  % desired 3D IV excitation
dIV_hi = f.iv;

f = matfile('./objSupportMask.mat');  % object support mask
m_hi = f.objSupportMask;

figure,
subplot(131), im(dIV_hi)
subplot(132), im(m_hi)
subplot(133), im(b0Map_hi)

% FOV and offset of `b0Map_hi`, `dIV_hi`, and `m_hi`.
[fov, offset] = deal([24, 24, 24], [0, 0, 0]);

doOV = true;
doSave = false;

%%
[pIV, pOV] = designScript(dIV_hi, m_hi, b0Map_hi, fov, offset, doOV, doSave);

% The outputs are structures with fields: `RF`, and `GR`
disp('pIV:')
disp(pIV)

disp('pOV:')
disp(pOV)
end
