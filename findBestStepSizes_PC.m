function [optimalSigma, optimalTau] = findBestStepSizes_PC(minSigma,... 
  maxSigma, minTau, maxTau, nrmK, sinogram,...
  nDetectors, detSize, thetas, translations, nCols, nRows, pixSize, ...
  varargin )

  defaultLevel=0;
  p = inputParser;
  p.addOptional( 'level', defaultLevel );
  p.parse( varargin{:} );
  level = p.Results.level;

  levelThresh = 3;
  nSteps = 10;

  sigmas = logspace(log10(minSigma),log10(maxSigma),nSteps);
  taus = logspace(log10(minTau),log10(maxTau),nSteps);

  nSigmas = numel(sigmas);
  nTaus = numel(taus);
  mincosts = 9999 * ones(nSigmas,nTaus);
  disp(['Working on level ', num2str(level)]);
  for m = 1:nSigmas
    disp(['findBestStepSizes: Working on ', num2str(m), ...
      ' of ', num2str(nSigmas)]);

    thisSigma = sigmas(m);

    parfor n = 1:nTaus
      % check if sigma*tau is within bound
      if thisSigma*taus(n) > 1 / (nrmK*nrmK), continue; end;

      [~, costs] = ctCorrectForTranslation_PC( sinogram, nDetectors, ...
        detSize, thetas, translations, nCols, nRows, pixSize, ...
        thisSigma, taus(n) );

      mincosts(m,n) = min( costs(:) );
    end
  end

  [minTaus, rowIndexes] = min(mincosts);
  [~, colIndex] = min(minTaus);

  optimalTau = taus(colIndex);
  rowSigma = rowIndexes(colIndex);
  optimalSigma = sigmas(rowSigma);

  level = level+1;

  if level <= levelThresh
    minSigma = sigmas(max(rowSigma-1,1));
    maxSigma = sigmas(min(rowSigma+1,nSteps));
    minTau = taus(max(colIndex-1,1));
    maxTau = taus(min(colIndex+1,nSteps));

    [optimalSigma, optimalTau] = findBestStepSizes_PC(minSigma,... 
      maxSigma, minTau, maxTau, nrmK, sinogram,...
      nDetectors, detSize, thetas, translations, nCols, nRows, ...
      pixSize, level );
  end

end
