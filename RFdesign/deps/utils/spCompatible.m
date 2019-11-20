function x = spCompatible(x)
% This function ensures a variable is compatible with matlab sparse, which can
% otherwise cause runtime error.

if issparse(x)||isa(x,'double')||islogical(x), return;
else                                         , x = double(x);
end

end

