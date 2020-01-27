%normalized inverse FFT
function [output] = Nifft2(matrix)
const = sqrt(size(matrix,1)*size(matrix,2));
output = ifft2(matrix)*const;
return;
