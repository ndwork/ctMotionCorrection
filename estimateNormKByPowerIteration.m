function [nrm,lambdaVals] = estimateNormKByPowerIteration( ...
  applyE,applyETrans,applyD1,applyD1Trans,applyD2,applyD2Trans, ...
  maxIters,x0)
  % x0 is an initial vector to get started with.

  if nargin < 8
    maxIters = 500;
  end

  stopTol = 1e-10;

  x = x0;
  lambdaPrev = 0;
  lambdaVals = zeros(maxIters,1);

  for iter = 1:maxIters
    disp(['Power Iteration: ', num2str(iter)]);

    D1x = applyD1(x);
    D2x = applyD2(x);
    D1TD1x = applyD1Trans(D1x);
    D2TD2x = applyD2Trans(D2x);
    Ex = applyE(x);
    ETEx = applyETrans(Ex);

    ATAx = ETEx + D1TD1x + D2TD2x;

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