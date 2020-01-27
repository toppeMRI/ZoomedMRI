% normalized FFT
function [output] = Nfft2(matrix)
const = sqrt(size(matrix,1)*size(matrix,2));
output = fft2(matrix)/const;
return;

