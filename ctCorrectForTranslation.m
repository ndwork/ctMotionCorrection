function [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize, varargin )
  % Optional Variables:
  % method:
	%   'GD' for Gradient Descent
  %   'LADMM' for Linearized ADMM
  %   'PC' for Pock Chambolle

  defaultMethod = 'PC';
  defaultR = [];
  p = inputParser;
  p.addOptional( 'method', defaultMethod );
  p.addOptional( 'radonMatrix', defaultR );
  p.parse( varargin{:} );
  method = p.Results.method;
  radonMatrix = p.Results.radonMatrix;

  switch method

    case 'GD'     % Gradient Descent
      costs = [];
      recon = ctCorrectForTranslation_GD( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, 'radonMatrix', radonMatrix );

    case 'LADMM'  % Linearized ADMM
      [recon,costs] = ctCorrectForTranslation_LADMM( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, 'radonMatrix', radonMatrix );

    case 'PC'     % Pock-Chambolle
      [recon,costs] = ctCorrectForTranslation_PC( sinogram, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, 'radonMatrix', radonMatrix );

  end

end

