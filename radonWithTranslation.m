
function sinogram = radonWithTranslation( img, pixSize, nDetectors, ...
  detSize, thetas, translations, verbose )
  % img:  2D array - will take the Radon transform of this image
  % pixSize: horizontal and vertical size of pixel (assumed square)
  % nDetectors: the number of detectors
  % thetas: an N element 1D array, each element is the angle that
  %   corresponds to row radon domain
  % translations is a Nx2 array specifying [y,x] meters translated
  %   per radian of rotation

  if nargin < 7, verbose = 0; end;
  
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
%   radiusImg = img .* radiusMask;
radiusImg = img;

  if mod( Nx, 2 )==0
    locs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * pixSize;
  else
    locs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * pixSize;
  end

  sinogram = zeros( nTheta, nDetectors );
  parfor th=1:numel(thetas)
    theta = thetas_deg(th);
    thisTrans_m = translations(th,:);
    thisTrans_pix = thisTrans_m / pixSize;
    translated = translateImg( radiusImg, thisTrans_pix );
    rotImg = imrotate( translated, theta, 'bilinear', 'crop' );
    sumResult = sum( rotImg, 1 ) * pixSize;

    interped = interp1( locs, sumResult, dLocs,'linear',0 );

    sinogram(th,:) = interped;

    if verbose && mod(th,10)==0
      disp(['radonWithTranslation Theta: ', num2str(th), ' of ', ...
        num2str(numel(thetas)) ]);
    end
  end

end

