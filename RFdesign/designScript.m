function [pIV, pOV] = designScript(dIV, m, b0Map, fov, offset, doOV, doSave)
% function [pIV, pOV] = designScript(dIV, m, b0Map, fov_ofst_cm, doSave)
%INPUTS
% - dIV (nx, ny, nz), hi-res desired IV excitation pattern
% - m (nx, ny, nz), hi-res object support
% - b0Map (nx, ny, nz), hi-res b0Map
% - fov: "cm", (1,3)
% - offset: "cm", (1,3)
% - doOV: [t/F], design OV pulse or not
% - doSave: [t/F], save design information or not
%OUTPUTS
% - pIV, structure, raw IV pulse:
%   .RF, (nT_IV, 1) pulse RF (complex)
%   .GR, (nT_IV, 3) pulse gradient
% - pOV, structure, raw OV pulse:
%   .RF, (nT_OV, 1) pulse RF (complex)
%   .GR, (nT_OV, 3) pulse gradient
%
% Pulse design script for zMRI

doCim = true;
imSize_lo = [32, 32, 20]; % SET MANUALLY, limited by RAM

if ~exist('offset', 'var'), offset = [0, 0, 0]; end
if ~exist('doOV', 'var'), doOV = false; end
if ~exist('doSave', 'var'), doSave = false; end

fn_hi2lo = @(X)cast(imresize3(double(X), imSize_lo, 'linear'), 'like', X);
fn_Mxy   = @(M3d) M3d(:,:,:,1) + 1i*M3d(:,:,:,2);
fn_cube2phm = @(cube,w)struct('sMap',cube.sMap, 'w',w, 'fov',cube.fov ...
                              , 'b0Map',cube.b0Map, 'ofst',cube.ofst);
fn_pulse = @(rf, g)struct('RF',rf, 'GR',g);

%% resize
dIV_lo = fn_hi2lo(dIV);
m_lo = fn_hi2lo(m);
b0Map_lo = fn_hi2lo(b0Map);

%% design pattern
dOV_lo = double(~imdilate(~~dIV_lo, ones(3,3,3)));

[mIV_lo, mOV_lo] = deal(m_lo&~~dIV_lo, m_lo&~~dOV_lo);

cube_lo = mCube(fov, imSize_lo, offset, 'm',true(imSize_lo), 'b0Map',b0Map_lo);

wIV_lo = zeros(imSize_lo);
[wIV_lo(mIV_lo), wIV_lo(mOV_lo)] = deal(2, 1);
phmIV_lo = fn_cube2phm(cube_lo, wIV_lo);

wOV_lo = zeros(imSize_lo);
[wOV_lo(mOV_lo), wOV_lo(mIV_lo)] = deal(2, 1);
phmOV_lo = fn_cube2phm(cube_lo, wOV_lo);

%% design the IV and OV pulses using Hao's method with size-scaled problem
[initMeth, isBalancedIV, isBalancedOV] = deal('xkt', true, true);

[methPara.npe, methPara.tol] = deal(90, 0.2);
[rfIV, gIV, ~] = ...
  pulse_st_sun(phmIV_lo, dIV_lo, isBalancedIV, initMeth, ...
               'methPara',methPara, 'ncycle',100, 'nIP',30);
drawnow; close(gcf); % TSP function will pop a irrelevant figure

pIV = fn_pulse(rfIV, gIV);

%% Display simulated single shot results
Miv = cube_lo.embed(cube_lo.applyPulse(pIV, doCim, false));
Miv_xy = fn_Mxy(Miv);

figure
subplot(121), im(Miv_xy); caxis([0,1]); colormap gray

%% Same procedure for OV
if doOV
  [methPara.npe, methPara.tol] = deal(110, 0.2);
  [rfOV, gOV, ~] = ...
    pulse_st_sun(phmOV_lo, dOV_lo, isBalancedOV, initMeth, ...
    'methPara',methPara, 'ncycle',100, 'nIP',30);
  drawnow; close(gcf); % TSP function will pop a irrelevant figure
  pOV = fn_pulse(rfOV, gOV);

  Mov = cube_lo.embed(cube_lo.applyPulse(pOV, doCim, false));
  Mov_xy = fn_Mxy(Mov);

  subplot(122), im(Mov_xy); caxis([0,1]); colormap gray
  drawnow
end

%% save design info
if ~doSave, return; end

system('rm designInfo.mat');

mfile = matfile('designInfo.mat');
mfile.b0Map = b0Map;
mfile.b0Map_lo = b0Map_lo;

mfile.wIV_lo = wIV_lo; % weighting for design
mfile.mIV_lo = mIV_lo; % support mask
mfile.dIV_lo = dIV_lo; % IV excitation pattern
mfile.Miv  = Miv; % raw IV pulse excitation profile, to be scaled

if doOV
  mfile.dOV_lo = dOV_lo; % OV excitation pattern
  mfile.Mov  = Mov; % raw OV pulse excitation profile, to be scaled
end

end

