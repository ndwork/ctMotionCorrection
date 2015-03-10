function recon = ctCorrectForRotAndTrans_GD( sinogram, nDetectors, ...
  detSize, thetas, rotations, translations_m, nCols, nRows, pixSize )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
    detSize, thetas );
  RT = transpose(R);

  translations_pix = translations_m / pixSize;

  applyE = @(u) radonWithRotAndTrans_GD( u, rotations, translations_pix, ...
    nDetectors, nCols, nRows, R );

  nThetas = numel( thetas );
  applyET = @(u) backprojectionWithRotAndTrans_GD( u, rotations, ...
    translations_pix, nCols, nDetectors, nThetas, RT );
  
  
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


function out = radonWithRotAndTrans_GD( u, rotations, translations_pix, ...
  nDetectors, nCols, nRows, R )

  img = reshape( u, [nRows nCols] );
  sinogram = RWithRotAndTrans( img, rotations, translations_pix, ...
    nDetectors, R );
  out = sinogram(:);
end

function out = backprojectionWithRotAndTrans_GD( u, rotations, ...
  translations_pix, nCols, nDetectors, nThetas, RT )

  sinogram = reshape( u, [nThetas nDetectors] );
  bp = RTWithRotAndTrans( sinogram, rotations, translations_pix, ...
    nCols, RT );

  out = bp(:);
end
