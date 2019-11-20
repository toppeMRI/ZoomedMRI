function isNonEmpty = nonEmptyVarChk(vn, allowEmpty)
% Dirty trick to test whether a variable exists and is not empty.
% One may forget sanity checking if an input is empty. An empty input will pass
% the exists check, but may crash a program.
%INPUTS
% - vn str, name of variable to be checked
% - allowEmpty [t/F], allow this var to be empty
%OUTPUTS
% - isvalid boolean, whether vn is a valid variable

if nargin == 0, test(); return; end

if ~exist('allowEmpty','var')||isempty(allowEmpty), allowEmpty=false; end

cmd = ['exist(''',vn,''',''var'')'];
if allowEmpty, isNonEmpty = evalin('caller', cmd);
else,          isNonEmpty = evalin('caller',[cmd,'&&~isempty(',vn,')']);
end

end

function test()

a = 1; %#ok<NASGU>
assert(nonEmptyVarChk('a') == true);
b = []; %#ok<NASGU>
assert(nonEmptyVarChk('b') == false);
assert(nonEmptyVarChk('b',true) == true);
disp('validVarChk.test() done');

end
