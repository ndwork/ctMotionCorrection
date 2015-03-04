function [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize, nrmK, varargin )
  % Optional Variables:
  % method:
	%   'GD' for Gradient Descent
  %   'LADMM' for Linearized ADMM
  %   'PC' for Pock Chambolle


  defaultMethod = 'PC';  
  p = inputParser;
  p.addOptional( 'method', defaultMethod );
  p.parse( varargin{:} );
  method = p.Results.method;

  switch method

    case 'GD'     % Gradient Descent
      costs = [];
      recon = ctCorrectForTranslation_GD( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize );

    case 'LADMM'  % Linearized ADMM
      [recon,costs] = ctCorrectForTranslation_LADMM( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, nrmK );

    case 'PC'     % Pock-Chambolle
      [recon,costs] = ctCorrectForTranslation_PC( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, nrmK );

  end

end

