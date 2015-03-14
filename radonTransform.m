
function sinogram = radonTransform( img, pixSize, nDetectors, ...
  detSize, thetas )
  % img:  2D array - will take the Radon transform of this image
  % pixSize: horizontal and vertical size of pixel (assumed square)
  % nDetectors: the number of detectors
  % thetas: an N element 1D array, each element is the angle that
  %   corresponds to row radon domain

  dOffset=0;  % center channel offset

  nTheta = numel(thetas);
  dLocs = ( [0:nDetectors-1] - 0.5*(nDetectors-1) ) * detSize - dOffset;
  thetas_deg = thetas * 180/pi;

  Ny = size( img, 1 );  halfY = (Ny+1)/2;
  Nx = size( img, 2 );  halfX = (Nx+1)/2;
  xs = ones(Ny,1) * (1:Nx);
  ys = (1:Ny)' * ones(1,Nx);
  xs = xs - halfX;
  ys = ys - halfY;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(Nx/2,Ny/2);
	radiusImg = img .* radiusMask;

  locs = ( (0:Nx-1) - 0.5*(Nx-1) ) * pixSize;

  M = makeLinearInterpMatrix( locs, dLocs );

  sinogram = zeros( nTheta, nDetectors );
  parfor th=1:numel(thetas)
    theta = thetas_deg(th);
    rotImg = imrotate( radiusImg, theta, 'bilinear', 'crop' );
    sumResult = sum( rotImg, 1 ) * pixSize;

    %extrapVal = 0;
    %interped = interp1( locs, sumResult, dLocs, 'linear', extrapVal );
    interped = M * sumResult';

    sinogram(th,:) = interped;
  end

end


