function [nrm,lambdaVals] = estimateNormByPowerIteration( ...
  applyA,applyATrans,x0,maxIters)
  % x0 is an initial vector to get started with.

  x = x0;
  
  if nargin < 4
    maxIters = 500;
  end

  stopTol = 1e-10;

  lambdaPrev = 0;
  lambdaVals = [];

  for iter = 1:maxIters

     ATransAx = applyATrans(applyA(x));
     
     lambda = norm(ATransAx(:));
     lambda = norm(applyA(x));
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