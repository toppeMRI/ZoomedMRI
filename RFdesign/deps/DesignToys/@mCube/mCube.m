% mCube Object:
%{
% Coded for excitation simulation in 3D
% 
% See Also: @mSpin
%
% Brief intro on methods:
%
%   Tianrui Luo, 2017
%}

classdef mCube < mSpin
  properties (SetAccess = immutable)
    fov; % (xyz,) cm, field of view
    loc; % (prod(dim), xyz) cm, loc of spins trade off RAM for easier coding
    ofst_cm; % (xyz, ) cm, offset of the fov;
  end
  properties (Dependent)
    ofst; % a.u. [-0.5, 0.5], ofst_cm./fov
    res; % (xyz,) cm, resolution
  end
  properties (SetAccess = public, GetAccess = public)
    dphs; % (dim), radians
  end
  
  methods (Access = public)
    function obj = mCube(fov, dim, ofst_cm, varargin)
      % varargin = {'w',w, 'b0Map',b0Map, 'sMap',sMap};
      obj@mSpin(dim, varargin{:});
      
      obj.fov = fov;
      mlocn = mask2mlocn(obj.m); % (nSpins, 1+xyz), 1+ for the column of mask
      obj.loc = bsxfun(@plus, bsxfun(@times, mlocn(:,2:end), fov), ofst_cm(:)');
      obj.ofst_cm = ofst_cm;
    end
    
    function [Mo, Mhst] = applyPulse(obj, pulse, doCim, doUpdate, dt)
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
      
      m = obj.m;
      Beff = obj.genBeff(pulse.GR, pulse.RF, obj.loc(m,:));
      
      [Mo, Mhst] = obj.applyBeff(Beff, doCim, doUpdate, dt);
    end
    
  end % of methods (Access = public)
  
  methods (Access = protected, Hidden = true)
  end % of methods (Access = private, hidden = true)
  
  methods % set & get methods
    function res = get.res(obj)
      res = obj.fov ./ obj.dim;
    end
    function ofst = get.ofst(obj)
      ofst = obj.ofst_cm./obj.fov;
    end
  end
  
  methods (Static = true)
  end
  
  methods (Access = private, Hidden = true) % obsolete or in-draft functions
    function output = imresize3(input,scale)
      % Modified from Hao's imresize3()
      if scale~=1
        dims = size(input);
        dims(1:2) = ceil(dims(1:2)*scale);
        output = cast(zeros(dims), 'like', input);
        for iz = 1:prod(dims(3:end)) % input might be more than 3d
          output(:,:,iz) = imresize(input(:,:,iz),scale, 'bilinear');
        end
      else
        output = input;
      end
    end
  end
end

% TODOs
% 1. modify to allow bx, by to be zero in applyPulse(), so AQ pulse can be
%    simed
% 2. should I enable randomess in SSFP_sim?
%

