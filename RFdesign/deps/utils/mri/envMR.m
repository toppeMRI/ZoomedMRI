function varargout = envMR(mod, varargin)
% contains environmental variables: gMax, sMax, dt, gam
%INPUTS:
% - mod str, case insensitive, {'set', 'get', 'get_s', 'reset'}
%   'set':   Set environmental persistent variables;
%   'get':   Get the values of persistent variables, return them separately;
%   'get_s': Get the values of queried variables, retuned as a structure;
%   'reset': Resetting the persistent variables to default values.
%NOTE:
% Many pulse design functions have separate environmental variable inputs that
% allows env vars to differ from the 'persistent' system env vars here.
% For instance, it can happen when one wants to empirically ensure the designed
% pulse does not violate any system env var restrictions.
% This heuristic of staying within the system limits is often due to the algs
% used are penalized methods rather than constraint methods.

if nargin == 0, test(); return; end

persistent mrEnv_p % _p: persistent, good to put these vars at the beginning.
fName_c = {'dt','gam','gMax','sMax'}; % {Sec, Hz/Gauss, Gauss/cm, Gauss/cm/Sec}
if isempty(mrEnv_p) || strcmpi(mod, 'reset')
  % GE MR750 setting used in fMRI lab UMich
  mrEnv_p = cell2struct({4e-6, 4257.6, 4, 12e3}, fName_c, 2);
  if strcmpi(mod, 'reset'), disp('...mrEnv_p in envMR reset');
  else, disp('...mrEnv_p in envMR initialized');
  end
end

% for cell input, e.g. envMR('get', {'dt', 'gam'}), envMR('set',{'dt'}).
if numel(varargin)==1 && iscell(varargin{1}), varargin = varargin{1}; end

switch lower(mod)
  case {'set', 'reset'}
    if ~strcmpi(mod,'reset') && ~isempty(varargin)
      mrEnv_p = attrParser(mrEnv_p, varargin, false);
    end
    varargout = {mrEnv_p};
  case 'get'
    varargout = getattrs(mrEnv_p, varargin, false);
    assert(nargout==numel(varargin));
  case 'get_s'
    if isempty(varargin), varargout = {mrEnv_p};
    else, varargout = {cell2struct(getattrs(mrEnv_p,varargin,0), varargin, 2)};
    end
  otherwise,  error('Unsupported set_get query');
end

end

function test()
prefix = mfilename('fullpath');
disp('------------------------');
disp([prefix, '.test()']);
env1_s = envMR('get_s'); %_s: struct, back up current values
fName_c = {'dt','gam','gMax','sMax'}; % {Sec, Hz/Gauss, Gauss/cm, Gauss/cm/Sec}
assert(isequal(fName_c, fieldnames(env1_s)'), 'env var name mismatch');

[dt, gam, gMax, sMax] = envMR('get', fName_c{:});
assert(isequal(getattrs(env1_s,fName_c),{dt,gam,gMax,sMax}), 'get test failed');

env0_s = envMR('reset');
env2_s = envMR('set', 'dt',exp(env0_s.dt)+exp(dt)); % env2.dt~= env[01].dt
assert(~isequal(env1_s,env2_s), 'set test failed');
assert(isequal(envMR('reset'), env0_s), 'reset test failed');

% restore from backup
envMR('set', env1_s);
disp([prefix, '.test() passed']);
end
