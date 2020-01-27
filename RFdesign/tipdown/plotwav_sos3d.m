
function plotwav(fname)
[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav(fname);
b1 = rho.*exp(i*theta);
plotrf(b1,gx,gy,gz);
