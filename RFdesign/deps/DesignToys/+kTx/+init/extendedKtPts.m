function [ko_pcm] = extendedKtPts(phm, d, npe, tol, gMax, sMax, dt)
% extended kt points
% Ref. Sun et al. DOI: 10.1109/TMI.2015.2478880
%
%EXTERNAL DEPENDENCIES
% tsp_ga: (MathWorks File Exchange)
%   goo.gl/kJQV6b
% minTimeGradient() (M. Lustig)
%   goo.gl/WfaJNU
%
%INPUTS
% - phm struct, with .sMap, .w (weight), .fov fields
% - d   (nx, ny, nz), desired excitation pattern
% - npe (1,) [70],  number of pe pts
% - tol (1,) [0.2], tolerance for min time gradient alg.
% - dt/gMax/sMax, hardware limits, see also envMR()
%OUTPUTS
% - ko_pcm (nk, nd), cycle/cm, extended kt_pts traj
%     not projected onto feasible k-space, be aware of this
%

[sMap, w, fov] = deal(phm.sMap, phm.w, phm.fov);
% parameter for cubic spline smoothing
if ~nonEmptyVarChk('npe'), npe = 70; end
if ~nonEmptyVarChk('tol'), tol = 0.2; end

[gMax0, sMax0, dt0] = envMR('get', 'gMax', 'sMax', 'dt');
if ~nonEmptyVarChk('gMax'), gMax = gMax0; end % Gauss/cm
if ~nonEmptyVarChk('sMax'), sMax = sMax0; end % Gauss/cm/Sec
if ~nonEmptyVarChk('dt'),   dt   = dt0;   end % Sec

%% units conversion
dt = dt*1e3; % mSec;
sMax = sMax*1e-3; % G/cm/Sec -> G/cm/mSec

%% find kt points, sub-routines locate in @PulseST/private
ncand = 10; % # candidates in KTpoints greedy search, empirically chosen (Sun)
k = getPePtsGreedy(d,npe,ncand,fov,w,sMap);

%% Use TSP and minTimeGrad
% TSP & minTimeGrad packed into @PulseST/private/TSPminTimeConnect.m,
% for re-usability;
ko_pcm = TSPminTimeGradConnect(k, tol, fov, dt, gMax, sMax);

end

function [ko_pcm] = TSPminTimeGradConnect(ki, tol, fov, dt, gMax, sMax)
% Need tsp_ga alg from mathwork
% Need minTimeGradient() from Lustig, M.
% knxyz(end,:) == koxyz(end,:), this method keeps end point of the input
% unchanged, i.e. it finds a shortest path getting throught all points and
% ends at the give end point
% Input
%  ki, (npts, xyz), cycle/FOV
% Output
%  ko, (nSteps, xyz), cycle/FOV

%% Re-Connect those points using tsp algorithm
% find the shortest path visiting order (original KT-points);
% from tsp_ga demo

npts = size(ki,1);
nIter      = 5e3;
popSize    = 20;
showProg   = 0;
showResult = 1;
a = meshgrid(1:npts);
dmat = reshape(sqrt(sum((ki(a,:)-ki(a',:)).^2,2)),npts,npts);

[depName, depURL] = deal('tsp_ga', 'goo.gl/kJQV6b');
if isempty(which(depName)), error(['Missing',depName, ', check',depURL]); end
optRoute = tsp_ga(ki,dmat,popSize,nIter,showProg,showResult);

kstart = find(optRoute == npts);
optRoute = [optRoute(kstart+1:end), optRoute(1:kstart)];
kn = ki(optRoute,:);

knorm = sum(kn.^2,2);
non0_ind = find(knorm~=0);
% below shortens the length of last kTrajectory segment
if knorm(non0_ind(end)) > knorm(non0_ind(1))
  kn((non0_ind(1):non0_ind(end)),:) = flipud(kn((non0_ind(1):non0_ind(end)),:));
end
% next line is just a left over from Hao, doesn't seems to matter
kn(1:max(non0_ind(1)-1,1), :) = [];

%% smooth the reconnected kt points
w = ones(size(kn,1),1);
w([1 end]) = 100; %smooth the kt points
kx = kn(:,1);    ky = kn(:,2);    kz = kn(:,3);
[~, kx_sp] = spaps(1:size(kn,1), kx.', tol, w, 3);
[~, ky_sp] = spaps(1:size(kn,1), ky.', tol, w, 3);
[~, kz_sp] = spaps(1:size(kn,1), kz.', tol, w, 3);
kn = [kx_sp; ky_sp; kz_sp].';

k_pcm = bsxfun(@rdivide, kn, fov);

%% Find the shortest time gradient given the path
% Rotationally invariant solution
rv = 0;
% minTimeGradient(): optimal control can have small violation of constraints.
[depName, depURL] = deal('minTimeGradient', 'goo.gl/WfaJNU');
if isempty(which(depName)), error(['Missing',depName, ', check',depURL]); end
[~, ~, ~, ~, K_rv] = minTimeGradient(k_pcm, rv,0,0,gMax,sMax, dt);

ko_pcm = K_rv;
end
