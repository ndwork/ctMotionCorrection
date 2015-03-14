
function recon = ctCorrectForTranslation_GD( sinogram, nDetectors, ...
  detSize, thetas, translations_m, nCols, nRows, pixSize, varargin )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  defaultR = [];
  p = inputParser;
  p.addOptional( 'radonMatrix', defaultR );
  p.parse( varargin{:} );
  R = p.Results.radonMatrix;

  if numel(R)==0
    R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
      detSize, thetas);
  end
  RT = transpose(R);

  translations_pix = translations_m / pixSize;

  applyE = @(u) radonWithTranslation_GD( u, translations_pix, ...
    nDetectors, nCols, nRows, RT );

  nThetas = numel( thetas );
  applyET = @(u) backprojectionWithTranslation_GD( u, translations_pix, ...
    nCols, nDetectors, nThetas, RT );
  
  
  function out= applyForGD( in, type )
    switch type
      case 'transp'
        out = applyET( in );
      case 'notransp'
        out = applyE( in );
    end
  end

  maxIter = 1000;
  lsqrOut = lsqr( @applyForGD, sinogram(:), 1e-06, maxIter );
  recon = reshape( lsqrOut, [nRows nCols] );

end


function out = radonWithTranslation_GD( u, translations_pix, ...
  nDetectors, nCols, nRows, RT )

  img = reshape( u, [nRows nCols] );
  sinogram = RWithTranslation( img, translations_pix, nDetectors, RT );
  out = sinogram(:);
end

function out = backprojectionWithTranslation_GD( u, translations_pix, ...
  nCols, nDetectors, nThetas, RT )

  sinogram = reshape( u, [nThetas nDetectors] );
  bp = RTWithTranslation( sinogram, translations_pix, nCols, RT );

  out = bp(:);
end


