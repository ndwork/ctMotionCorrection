function [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize, nrmK, ...
  sigma, tau, varargin )

  defaultMethod = 'linearizedADMM';  
  p = inputParser;
  p.addOptional( 'method', defaultMethod );
  p.parse( varargin{:} );
  method = p.Results.method;

  switch method

    case 'linearizedADMM'
      [recon,costs] = ctCorrectForTranslation_LADMM( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, nrmK );

    case 'PC'
      [recon,costs] = ctCorrectForTranslation_PC( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, nrmK, sigma, tau );

  end

end

