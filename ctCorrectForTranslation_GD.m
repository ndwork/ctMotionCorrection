
function recon = ctCorrectForTranslation_GD( sinogram, nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  if nargin < 10, lambda=1; end;

  applyE = @(u) radonWithTranslation_GD( u, nCols, nRows, pixSize, ...
    nDetectors, detSize, thetas, translations );

  cx = 0;  cy = 0;
  applyET = @(u) backprojectionWithTranslation_GD( u, nDetectors, ...
    thetas, detSize, cx, cy, nCols, nRows, pixSize, translations );
  
  
  function out= applyForGD( in, type )
    switch type
      case 'transp'
        out = applyET( in );
      case 'notransp'
        out = applyE( in );
    end
  end
  
  lsqrOut = lsqr( @applyForGD, sinogram(:) );
  recon = reshape( lsqrOut, [nRows nCols] );

end


function out = radonWithTranslation_GD( u, nCols, nRows, pixSize, ...
  nDetectors, detSize, thetas, translations )

  img = reshape( u, [nRows nCols] );

  sinogram = radonWithTranslation( img, pixSize, nDetectors, ...
    detSize, thetas, translations );
  
  out = sinogram(:);
end

function out = backprojectionWithTranslation_GD( u, nDetectors, thetas, ...
  detSize, cx, cy, nCols, nRows, pixSize, translations )
  
  nThetas = numel( thetas );
  sinogram = reshape( u, [nThetas nDetectors] );

  %bp = radonWithTranslationAdjoint( sinogram, thetas, ...
  %  detSize, cx, cy, nCols, nRows, pixSize, translations );
  bp = backprojectionWithTranslation( sinogram, thetas, ...
    detSize, cx, cy, nCols, nRows, pixSize, translations );

  out = bp(:);
end


