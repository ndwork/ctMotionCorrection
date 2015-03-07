
function [recon,costs] = ctCorrectForRotAndTrans( sinogram, nDetectors, ...
  detSize, thetas, rotations, translations, nCols, nRows, pixSize, ...
  varargin )
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
      recon = ctCorrectForRotAndTrans_GD( sinogram, ...
        nDetectors, detSize, thetas, rotations, translations, ...
        nCols, nRows, pixSize );

    case 'LADMM'  % Linearized ADMM
      [recon,costs] = ctCorrectForRotAndTrans_LADMM( sinogram, ...
        nDetectors, detSize, thetas, rotations, translations, ...
        nCols, nRows, pixSize );

    case 'PC'     % Pock-Chambolle
      [recon,costs] = ctCorrectForRotAndTrans_PC( sinogram, ...
        nDetectors, detSize, thetas, rotations, translations, ...
        nCols, nRows, pixSize );

  end

end

