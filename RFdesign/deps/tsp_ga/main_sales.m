n = 50;
xyz = 10*rand(n,3);
popSize = 60;
numIter = 1e4;
showProg = 1;
showResult = 1;
a = meshgrid(1:n);
dmat = reshape(sqrt(sum((xyz(a,:)-xyz(a',:)).^2,2)),n,n);
[optRoute,minDist] = tsp_ga(xyz,dmat,popSize,numIter,showProg,showResult);