% load the converge plot of all methods
if 0
   load weighting;
   beta = 8; 
   exetimeAll = [];
   nrmseAll = [];
   costAll = [];
   for i = 1:4
      fname = ['convergeSpeed', num2str(i)];
      load(fname);
      for iter = 1:length(k_allIter)
         k = k_allIter{iter}; 
         b = b_allIter{iter};
         A = formA(k, sens, mask, fieldmap, tt, xyzRange);
         betavec = sqrt(beta)*ones(size(A,2),1);
         W = diag_sp(weightIm(mask));
         b = qpwls_pcg1_hao(zeros(size(betavec)), A, W, d(mask), diag_sp(betavec),'niter',100);
         [e, As] = eval_func(A, b, d, sens, mask);
         nrmseFull(1,iter) = norm(e)/sqrt(numel(e))/max(abs(d(mask)));
      end
      exetimeAll(i,:) = exetime;
      nrmseAll(i,:) = nrmseFull;
      costAll(i,:) = cost;
   end
%    fname = ['convergeSpeed', num2str(4)];
%    load(fname);
   save plotConvergeData exetimeAll nrmseAll costAll exetime nrmse cost
else
   load plotConvergeData;
   
   %%%%% include the error after B-spline fitting
   load plotConvergeData1;
   exetimeAll = exetimeAll1; 
   nrmseAll = nrmseAll1; 
   costAll = costAll1; 
end

niter = 10;
figure,plot(exetimeAll(4,1:niter),nrmseAll(4,1:niter),'bv-'); 
hold on; plot(exetimeAll(1,:),nrmseAll(1,:),'+-','Color',[0,0.5,0]); 
hold on; plot(exetimeAll(2,:),nrmseAll(2,:),'ro-'); 
hold on; plot(exetimeAll(3,:), nrmseAll(3,:),'kx-');
xlim([0 32]); ylim([0.1, 0.135]);  
legend('fmincon', 'projected GD', 'projected LM', 'interior point');
xlabel('time (sec)'); ylabel('NRMSE');

figure,plot(exetimeAll(4,1:niter),costAll(4,1:niter),'bv-'); 
hold on; plot(exetimeAll(1,:),costAll(1,:),'+-','Color',[0,0.5,0]); 
hold on; plot(exetimeAll(2,:),costAll(2,:),'ro-'); 
hold on; plot(exetimeAll(3,:), costAll(3,:),'kx-');
xlim([0 32]); 
legend('fmincon', 'projected GD', 'projected LM', 'interior point');
xlabel('time (sec)'); ylabel('Cost function value (a.u.)');