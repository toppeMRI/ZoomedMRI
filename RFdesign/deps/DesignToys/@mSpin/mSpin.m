% mSpin Object:
% objects of spin isochromats
% This is a metaClass for mCube and mVoxels
%
% Brief intro on methods:
%
%   Tianrui Luo, 2017
%}

classdef mSpin < matlab.mixin.SetGet
  properties (SetAccess = immutable)
    dim;
  end
  properties (SetAccess = public, GetAccess = public)
    spins = []; % (nSpins, xyz); spin vector, rotating frame
    
    gam = [];   % Hz/G
    T1  = 1.47;    % Sec, dflt for grey matter
    T2  = 0.07;  % Sec, dflt for grey matter

    m     = []; % (nSpins, 1) a.u., simulation mask,
    b0Map = []; % (nSpins, 1) Hz
    sMap  = []; % (a,b,c, ncoil) a.u., transmit coil sMap, a*b*c == nSpins
    
    args  = []; % arbitrary contents user want to attach
  end
  
  methods (Sealed = true)
    function obj = mSpin(dim, varargin)
      % The constructor of a class will automatically create a default obj.
      st = cutattrs(attrParser(obj, varargin, false), properties('mSpin'));
      depProp_c = intersect(attrs(obj,'Dependent'), attrs(st));
      st = rmattrs(st, depProp_c); % rm dep props, which cannot be set
      for fName = fieldnames(st)', obj.(fName{1}) = st.(fName{1}); end
      
      if     isscalar(dim), dim = [dim, 1, 1];
      elseif numel(dim)==2, dim = [dim, 1];
      end
      obj.dim = dim; % Ensure required input is set correctly
      
      nSpins = prod(dim);
      if isempty(obj.gam),   obj.gam = envMR('get','gam'); end
      if isempty(obj.m),     obj.m = true(dim); end
      if isempty(obj.spins), obj.spins = [zeros(nSpins,2), ones(nSpins,1)]; end
    end
    
    function [Mo, Mhst] = blochSim(obj, beff, dt, doCim)
      %INPUTS:
      % - beff  (nSpins(masked), xyz) Gauss
      % - dt    (1,) Sec
      % - doCim (T/f), use C blochCim
      if ~nonEmptyVarChk('dt'), dt = envMR('get', 'dt'); end % Sec
      if ~nonEmptyVarChk('doCim'), doCim = true; end

      m_ = obj.m;
      if all(m_(:)), m_ = ':'; end % indexing speed trick
      [T1_, T2_, gam_] = deal(obj.T1, obj.T2, obj.gam);

      Mi = obj.spins(m_,:);
      [Mo, Mhst] = blochSim(Mi, beff, T1_, T2_, dt, gam_, doCim);
    end
    
    function Beff = genBeff(obj, GR, RF, loc)
      %Inputs:
      % - GR  (nSteps, xyz), Gauss/cm
      % - RF  (nSpins, ncoil, nSteps) or (nSteps,1), Gauss
      % - loc (nSpins[masked], xyz(, nStep)), cm
      %Outputs:
      % - Beff (nSpins[masked] xyz(, nStep)), effective field in rotating frame
      m_ = obj.m;
      if size(loc,1) == numel(m_), loc = loc(m_,:,:); end
      if iscolumn(RF), RF = reshape(RF, 1,1,[]); end
      if size(RF,1) == numel(m_), RF = RF(m_,:,:); end
      
      if ~isempty(obj.sMap)
        s = reshape(obj.sMap, [], size(obj.sMap, 4)); % -> ([], nc)
        if size(RF, 2) == 1, s = sum(s,2); end
        RF = bsxfun(@times, s(m_,:), RF);
      end
      RF = sum(RF, 2); % Gauss
      if size(RF,1) == 1, RF = repmat(RF, nnz(m_),1,1); end
      
      bz = reshape(loc*GR', nnz(m_), 1, []); % (nSpins[masked], 1, nSteps)
      % count in the B0 fieldmap effect
      dbz = 0;
      if ~isempty(obj.b0Map), dbz = obj.b0Map(m_)./obj.gam; end % Gauss
      bz = bsxfun(@plus, bz, dbz);
      
      Beff = [real(RF), imag(RF), bz]; % (nSpins[masked], xyz, nSteps)
    end
    
    function [Mo, Mhst] = applyBeff(obj, Beff, doCim, doUpdate, dt)
      % just a wrap of blochSim
      %INPUTS:
      % - pulse: container with following fields
      %   .RF (nSteps, nc), Gauss, complex, r(RF) = bx, i(RF) = by;
      %   .GR (nSteps, xyz), Gauss/cm
      % - dt (1,), Sec
      %OUTPUTS:
      % - Mhst (nspins(masked), xyz, nsteps), evolving history
      if ~nonEmptyVarChk('dt'),       dt = envMR('get', 'dt'); end % Sec
      if ~nonEmptyVarChk('doCim'),    doCim = true; end
      if ~nonEmptyVarChk('doUpdate'), doUpdate = true; end
      
      [Mo, Mhst] = obj.blochSim(Beff, dt, doCim);
      
      if doUpdate, obj.spins(obj.m,:) = Mo; end
    end
    
    function Mo = freePrec(obj, dur, doUpdate)
      % ROTATING FRAME
      % dur: Sec
      %INPUTS:
      % - dur (1,) Sec, duration of free precision
      % - doUpdate
      %OUTPUTS:
      % - spinVec (nSpins(masked), xyz)
      if ~nonEmptyVarChk('doUpdate'), doUpdate = true; end
      
      [Mo, b0Map_, T1_, T2_]= obj.mask(obj.spins, obj.b0Map, obj.T1, obj.T2);
      Mo = freePrec(Mo, dur, T1_, T2_, b0Map_);
      
      if doUpdate, obj.spins(obj.m,:) = Mo; end
    end
    
  end
  
  methods (Static = true)
  end
  
  methods % Utilities
    function varargout = mask(obj, varargin)
      m_ = obj.m;
      if all(m_(:)), m_ = ':'; end % indexing speed trick
      varargout = cell(1,numel(varargin));
      for ic = 1:numel(varargin)
        if isscalar(varargin{ic}), varargout{ic} = varargin{ic};
        else,                      varargout{ic} = varargin{ic}(m_,:);
        end
      end
    end

    function varargout = embed(obj, varargin)
      [m_, dim_] = deal(obj.m, obj.dim);
      if all(m_(:))
        varargout = cellfun(@(x)reshape(x,[dim_,size(x,2)]), varargin, 'Uni',0);
        return;
      end
      ne = prod(dim_); % #elements
      varargout = cell(1,numel(varargin));
      for ic = 1:numel(varargin)
        [~, nd] = size(varargin{ic});
        varargout{ic}       = zeros([ne, nd]);
        varargout{ic}(m_,:) = varargin{ic};
        varargout{ic}       = reshape(varargout{ic}, [dim_, nd]);
      end
    end
    
  end
  
  methods % set and get, sealed if the property cannot be redefined
    function set.spins(obj, spins)
      % Accepting spins organized in [nspins-xyz].
      eID = 'mSpin:set.spins:';
      assert(isreal(spins) && ~any(isnan(spins(:))), [eID,'get real']);
      if ~isempty(obj.spins)
        assert(isequal(size(spins),size(obj.spins)), [eID,'size mismatch']);
      end
      % self-multiplication is faster than power
      spins_n = sum(spins.*spins, 2);
      gt1 = spins_n>1;
      if any(gt1)
        warning('capped to unit-norm');
        spins(gt1) = spins(gt1)./sqrt(spins_n(gt1));
      end
      obj.spins = spins;
    end

  end

  methods (Sealed = true, Access = protected)
  end
  
  methods (Hidden = true, Access = private) % obsolete or in-draft fns
    function angleH =angHisFM(obj, dt, beff, spinTraj, phaseSq, iswrapped)
      % Input:
      %  dt, mSec
      %  beff, Tesla
      %    [nSpins(masked) xyz nSteps], nSpins > 1
      %    [nSteps xyz],        nSpins ==1
      %  spinTraj,  [nSpins(masked) xyz nSteps]
      %  phaseSq,  [nSteps 1], rad, describing phase of new x'-axis
      %  iswrapped, bool, wether phaseSq is wraped
      % Outputs
      %  angleH, angle between a spin and its corresponding beffFM
      %    [nSpins, nSteps]
      
      if size(beff,3) == 1 % assuming there is no one-step simu
        beff = reshape(beff', 1,3,[]); % should be faster than permute
        disp('beff is non spin-specific, assuming general');
      end
      
      dbz = obj.dbzFM(dt, phaseSq, iswrapped);
      bz = bsxfun(@plus, beff(:,3,:), dbz);
      beff = [beff(:,1:2,:), bz];
      
      angleH = obj.angHis(beff, spinTraj);
    end
    
    function angleH = angHis(beff, spinTraj)
      % Inputs
      %  beff, Tesla
      %    [nSpins(masked) xyz nSteps], nSpins > 1
      %    [nSteps xyz],        nSpins ==1
      %  spinTraj,  [nSpins(maksed) xyz nSteps]
      % Outputs
      %  angleH, rad, angle between a spin and its corresponding beff
      %    [nSpins, nSteps]
      if size(beff,3) == 1 % assuming there is no one-step simu
        beff = reshape(beff', 1,3,[]); % should be faster than permute
        disp('beff is non spin-specific, assuming general');
      end
      innerProd = sum(bsxfun(@times, beff, spinTraj), 2);
      beffNorm = sqrt(sum(beff.*beff, 2));
      spinTrajNorm = sqrt(sum(spinTraj.*spinTraj, 2));
      
      denom = bsxfun(@times, beffNorm, spinTrajNorm);
      
      angleH = acosd(innerProd ./ denom);
      if any(isnan(angleH(:)))
        warning('NaN angles are set 0');
      end
      angleH(isnan(angleH)) = 0;
      angleH = reshape(angleH, size(angleH,1), []);
    end

    function vecXYpTraj = yawFrame(phaseSq, vecXYTraj)
      % view vec3dTraj from frame yawed by phaseSq
      % Input:
      %  phaseSq, [nSteps 1], rad, describing phase of new x'-axis
      %  vec3dTraj, [nVector xy nsteps], vector Traj history
      % Output:
      %  vec3dTrajYawed, yawed vec3dTraj
      xTraj = vecXYTraj(:,:,1);
      yTraj = vecXYTraj(:,:,2);
      
      xyTraj = xTraj + 1i*yTraj;
      phaseSq = exp(1i*reshape(phaseSq,1,1,[]));
      % p for prime
      xypTraj = bsxfun(@times, xyTraj, -phaseSq); % minus phase change
      xpTraj = real(xypTraj);
      ypTraj = imag(xypTraj);
      vecXYpTraj = [xpTraj, ypTraj];
    end
    
    function beffFM = beffR2FM(obj, dt, beff, phaseSq, iswrapped)
      % beff rotation frame to frequency modulated frame
      % Inputs:
      %  dt, mSec
      %  beff, Tesla
      %   [nSpins(masked) xyz nSteps]
      %   [nSteps xyz]
      %  phaseSeq, [nSteps 1], rad, describing phase of new x'-axis
      %  iswrapped, bool, wether phaseSeq is wraped
      % Outputs:
      %  beffFM, Tesla, Frequency Modulated beff
      %   [nSpins(masked) xyz nSteps]
      %   [nSteps xyz]
      if ~nonEmptyVarChk('iswrapped'), iswrapped = true; end
      if size(beff,3) == 1 % assuming there is no one-step simu
        beff = reshape(beff', 1,3,[]); % should be faster than permute
        disp('beff is non spin-specific, assuming general');
      end
      
      dbz = obj.dbzFM(dt, phaseSq, iswrapped);
      bzFM = bsxfun(@plus, beff(:,3,:), dbz);
      beffFM = [yawFrame(phaseSq, beff(:,1:2,:)), bzFM];
      
      if size(beffFM,1) == 1 % assuming there is no one-step simu
        beffFM = reshape(beffFM, 3,[])'; % should be faster than permute
        disp('beff is non spin-specific');
      end
    end
    
    function dbz = dbzFM(obj, dt, phaseSq, iswrapped)
      % Inputs:
      %  dt, mSec
      %  phaseSq, [nSteps 1], rad, describing phase of new x'-axis
      %  iswrapped, bool, wether phaseSq is wraped
      % Outputs:
      %  dbz, [1 1 nSteps], Telsla
      if ~nonEmptyVarChk('iswrapped'), iswrapped = true; end
      gambar = 2*pi*obj.gam*1e4; % -> 2pi * Hz/T
      if iswrapped == true, phaseSq = unwrap(phaseSq); end
      
      dw = diff([0;phaseSq]) /dt *1e3;  % -> 2pi * Hz
      dbz = reshape(dw/gambar, 1,1,[]); % -> T
    end
    
  end
  
end


