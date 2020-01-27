function C = Cfd1_zeros(N)
% generate 1d finite differencing object with 0 boundary conditions. Notice that the
% last number of k is 0. 

forwc = @(arg, x) [x(2:end); 0] - x;
backc = @(arg, x) [0; x(1:end-1)] - x;
C = fatrix2('forw', forwc, 'back', backc, 'idim', N, 'odim', N);
