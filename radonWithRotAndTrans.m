
function sinogram = radonWithRotAndTrans( img, pixSize, nDetectors, ...
  detSize, thetas, rotations, translations_m )
  % img:  2D array - will take the Radon transform of this image
  % pixSize: horizontal and vertical size of pixel (assumed square)
  % nDetectors: the number of detectors
  % thetas: an N element 1D array, each element is the angle that
  %   corresponds to row radon domain
  % rotations: an N element 1D array, each element is the angle that
  %   the object is rotated
  % translations is a Nx2 array specifying [y,x] meters translated
  %   per radian of rotation

  dOffset=0;  % center channel offset

  nTheta = numel(thetas);
  dLocs = ( [0:nDetectors-1] - 0.5*(nDetectors-1) ) * detSize - dOffset;
  thetas_deg = thetas * 180/pi;

  Ny = size( img, 1 );  halfY = Ny/2;
  Nx = size( img, 2 );  halfX = Nx/2;
  xs = ones(Ny,1) * (1:Nx);
  ys = (1:Ny)' * ones(1,Nx);
  xs = xs - halfX;
  ys = ys - halfY;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(Nx/2,Ny/2);

  if mod( Nx, 2 )==0
    locs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * pixSize;
  else
    locs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * pixSize;
  end

  translations_pix = translations_m / pixSize;

  M = makeLinearInterpMatrix( locs, dLocs );

  rotations_deg = rotations * 180/pi;
  
  sinogram = zeros( nTheta, nDetectors );
  parfor th=1:numel(thetas)
    thisRot = rotations_deg(th);
    rotated = imrotate( img, thisRot, 'bilinear', 'crop' );
    
    thisTrans_pix = translations_pix(th,:);
    translated = translateImg( rotated, thisTrans_pix );
    
    masked = translated .* radiusMask;
    
    theta = thetas_deg(th);
    rotImg = imrotate( masked, theta, 'bilinear', 'crop' );
    sumResult = sum( rotImg, 1 ) * pixSize;

    extrapVal = 0;
    %interped = interp1( locs, sumResult, dLocs, 'linear', extrapVal );
    interped = M * sumResult';

    sinogram(th,:) = interped;
  end

end

