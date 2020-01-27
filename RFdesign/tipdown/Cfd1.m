function C = Cfd1(N)
% generate 1d finite differencing object with circulant boundary conditions

forwc = @(arg, x) circshift(x, 1) - x;
backc = @(arg, x) circshift(x, -1) - x;
C = fatrix2('forw', forwc, 'back', backc, 'idim', N, 'odim', N);
