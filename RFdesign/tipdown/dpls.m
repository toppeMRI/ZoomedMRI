function b = dpls(A,d,beta,mask)
% direct method for penalized least square
% solve ||Ab - d||_2^2 + beta ||b||_2^2

Atild = [A; sqrt(beta)*eye(size(A,2))];
dtild = [d(mask); zeros(size(A,2),1)]; 
b = Atild\dtild;

end
