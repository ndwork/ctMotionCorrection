function [optimalLambda, optimalMu] = findGoodStepSizes_LADMM(minLambda,... 
  maxLambda, nrmK, sinogram, nDetectors, detSize, thetas, translations, ...
  nCols, nRows, pixSize, varargin )

  defaultLevel=0;
  p = inputParser;
  p.addOptional( 'level', defaultLevel );
  p.parse( varargin{:} );
  level = p.Results.level;

  levelThresh = 3;
  nSteps = 10;

  lambdas = logspace(log10(minLambda),log10(maxLambda),nSteps);

  nLambdas = numel(lambdas);
  mincosts = 99999 * ones(nLambdas,1);
  disp(['Working on level ', num2str(level)]);
  for i = 1:nLambdas
    disp(['findGoodStepSizes: Working on ', num2str(i), ...
      ' of ', num2str(nLambdas)]);

    thisLambda = lambdas(i);
    thisMu = thisLambda / (nrmK*nrmK) * 0.999;

    [~, costs] = ctCorrectForTranslation_LADMM( sinogram, nDetectors, ...
      detSize, thetas, translations, nCols, nRows, pixSize, ...
      thisLambda, thisMu );

    mincosts(i) = min( costs(:) );
  end

  [~, minIndx] = min(mincosts);

  optimalLambda = lambdas( minIndx );
  optimalMu = optimalLambda / (nrmK*nrmK);

  level = level+1;
  if level <= levelThresh
    minLambda = lambdas(max(minIndx-1,1));
    maxLambda = lambdas(min(minIndx+1,nSteps));

    [optimalLambda, optimalMu] = findGoodStepSizes_LADMM( minLambda,... 
      maxLambda, nrmK, sinogram, nDetectors, detSize, thetas, ...
      translations, nCols, nRows, pixSize, level );
  end

end
