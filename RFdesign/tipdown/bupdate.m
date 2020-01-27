function b = bupdate(b, As, W, d, mask, betavec, niter)

b = qpwls_pcg1_hao(b, As, W, d(mask), diag_sp(betavec),'niter',niter);

%[b,FVAL,EXITFLAG,OUTPUT] = quadprog(As, dkB, Call, 0.99*(bound-Call*coeff_y),[],[],[],[],[],quadOptions);

