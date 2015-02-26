function [optimalSigma, optimalTau] = findBestStepSizes(minSigma,... 
  maxSigma, minTau, maxTau, nrmK, sinogram,...
  nDetectors, detSize, thetas, translations, nCols, nRows, pixSize, ...
  level )

  if nargin < 16, level=0; end;

  sigmas = logspace(log10(minSigma),log10(maxSigma),5);
  taus = logspace(log10(minTau),log10(maxTau),5);
  
%   mincosts = 9999 * ones(numel(sigmas),numel(taus));
  mincosts = [9999 0 0];
  for m = 1:numel(sigmas)
    parfor n = 1:numel(taus)
      % check if sigma*tau is within bound
      if sigmas(m)*taus(n) > 1 / (nrmK*nrmK), 
        continue;
      end;
      
      [~, costs] = ctCorrectForTranslation( sinogram, nDetectors, ...
        detSize, thetas, translations, nCols, nRows, pixSize, nrmK, ...
        sigmas(m), taus(n) );

%       mincosts(m,n) = min(costs(:));
        mincosts(m,n) = min( costs(:) );
    end
  end

  [minTaus, rowIndexes] = min(mincosts);
  [~, colIndex] = min(minTaus);
  
  optimalTau = taus(colIndex);
  optimalSigma = sigmas(rowIndexes(colIndex));

  level = level+1;
  
  if level < levelThresh
    minSigma = optimalSigma / 10;
    maxSigma = opticalSigma * 10;
    minTau = optimalTau / 10;
    maxTau = optimalTau * 10;
    findBestStepSizes(minSigma,... 
      maxSigma, minTau, maxTau, nrmK, sinogram,...
      nDetectors, detSize, thetas, translations, nCols, nRows, pixSize, ...
      level );
  end
  
  
end