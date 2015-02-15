function [nrm,lambdaVals] = estimateNormByPowerIterationED( ...
  applyE,applyETrans,applyD1,applyD1Trans,applyD2,applyD2Trans,x0,maxIters)
  % x0 is an initial vector to get started with.

  if nargin < 8
    maxIters = 500;
  end

  stopTol = 1e-10;

  x = x0;
  lambdaPrev = 0;
  lambdaVals = zeros(maxIters,1);

  for iter = 1:maxIters

    D1x = applyD1(x);
    D2x = applyD2(x);
    D1Tx = applyD1Trans(x);
    D2Tx = applyD2Trans(x);
    Ex = applyE(x);
    ETx = applyETrans(x);

    ATAx = ETx*Ex + D1Tx*D1x + D2Tx*D2x;

    lambda = norm(ATAx(:));
    lambdaVals(iter) = lambda;

    if lambda == 0, break; end;

    x = ATAx / lambda;

    diff = abs(lambda - lambdaPrev);
    if diff < stopTol, break; end;

    lambdaPrev = lambda;
  end

  lambdaVals = lambdaVals(1:iter);
  nrm = sqrt(lambda);

end