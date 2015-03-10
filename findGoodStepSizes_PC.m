function [optimalSigma, optimalTau] = findGoodStepSizes_PC(minSigma,... 
  maxSigma, nrmK, sinogram, nDetectors, detSize, thetas, translations, ...
  nCols, nRows, pixSize, varargin )

  defaultLevel=0;
  p = inputParser;
  p.addOptional( 'level', defaultLevel );
  p.parse( varargin{:} );
  level = p.Results.level;

  levelThresh = 3;
  nSteps = 10;

  sigmas = logspace(log10(minSigma),log10(maxSigma),nSteps);

  nSigmas = numel(sigmas);
  mincosts = 9999 * ones(nSigmas,1);
  disp(['Working on level ', num2str(level)]);
  for i = 1:nSigmas
    disp(['findGoodStepSizes: Working on ', num2str(i), ...
      ' of ', num2str(nSigmas)]);

    thisSigma = sigmas(i);
    thisTau = 1 / ( nrmK*nrmK*thisSigma ) * 0.999;

    [~, costs] = ctCorrectForTranslation_PC( sinogram, nDetectors, ...
      detSize, thetas, translations, nCols, nRows, pixSize, ...
      thisSigma, thisTau );

    mincosts(i) = min( costs(:) );
  end

  [~, minIndx] = min(mincosts);

  optimalSigma = sigmas( minIndx );
  optimalTau = 1 / ( nrmK*nrmK*optimalSigma ) * 0.999;

  level = level+1;
  if level <= levelThresh
    minSigma = sigmas(max(minIndx-1,1));
    maxSigma = sigmas(min(minIndx+1,nSteps));

    [optimalSigma, optimalTau] = findGoodStepSizes_PC( minSigma,... 
      maxSigma, nrmK, sinogram, nDetectors, detSize, thetas, ...
      translations, nCols, nRows, pixSize, level );
  end

end
