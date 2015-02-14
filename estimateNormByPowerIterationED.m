function [nrm,lambdaVals] = estimateNormByPowerIterationED( ...
  applyE,applyETrans,applyD1,applyD1Trans,applyD2,applyD2Trans,x0,maxIters)
  % x0 is an initial vector to get started with.

  x = x0;
  
  if nargin < 4
    maxIters = 500;
  end

  stopTol = 1e-10;

  lambdaPrev = 0;
  lambdaVals = [];

  for iter = 1:maxIters
    
      D1x = applyD1(x);
      D2x = applyD2(x);
      D1xTrans = applyD1Trans(x);
      D2xTrans = applyD2Trans(x);
      Ex = applyE(x);
      ExTrans = applyETrans(x);
      
      ATransAx = (ExTrans*Ex + D1xTrans*D1x + D2xTrans*D2x) * x;

%      ATransAx = applyATrans(applyA(x));
     
     lambda = norm(ATransAx(:));
     
     if lambda == 0
         break
     end
     x = ATransAx/lambda;

     diff = abs(lambda - lambdaPrev);
     if diff < stopTol
         break
     end

     lambdaPrev = lambda;  
     lambdaVals = [lambdaVals,lambda];

  end

  nrm = sqrt(lambda);


end