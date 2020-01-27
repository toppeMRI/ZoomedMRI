function roi = autoROI(b0map,m)
% b0map = Hz
% m = magnitude image

m = abs(m);
roi = m > 0.1*max(m(:)); % & ~edge(b0map,'canny');

return;

