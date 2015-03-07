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

  %if nargin < 16, level=0; end;

  sigmas = logspace(log10(minSigma),log10(maxSigma),5);
  taus = logspace(log10(minTau),log10(maxTau),5);

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
  optimalSigma = sigmas(rowIndexes(colIndex));

  level = level+1;

  if level <= levelThresh
    minSigma = optimalSigma / 10;
    maxSigma = optimalSigma * 10;
    minTau = optimalTau / 10;
    maxTau = optimalTau * 10;
    [optimalSigma, optimalTau] = findBestStepSizes_PC(minSigma,... 
      maxSigma, minTau, maxTau, nrmK, sinogram,...
      nDetectors, detSize, thetas, translations, nCols, nRows, ...
      pixSize, level );
  end

end
