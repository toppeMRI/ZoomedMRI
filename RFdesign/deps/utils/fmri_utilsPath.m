function path = genericGotoPath()
% This function name is a fake, it will be called by the file name, as from
% matlab document

filename = mfilename;
fullpath = mfilename('fullpath');

path = fullpath;
path(end-length(filename):end) = [];

if nargout == 0
  cd(path);
end

end
