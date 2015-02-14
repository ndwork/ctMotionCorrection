
function sinogram = radonWithTranslation( img, delta, nDetectors, ...
  dSize, thetas, translations )
  % img:  2D array - will take the Radon transform of this image
  % delta: horizontal and vertical size of pixel (assumed square)
  % nDetectors: the number of detectors
  % thetas: an N element 1D array, each element is the angle that
  %   corresponds to row radon domain
  % translations is a Nx2 array specifying [y,x] meters translated
  %   per radian of rotation

  dOffset=0;  % center channel offset

  nTheta = numel(thetas);

  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  thetas_deg = thetas * 180/pi;

  Ny = size( img, 1 );  halfY = Ny/2;
  Nx = size( img, 2 );  halfX = Nx/2;
  xs = ones(Ny,1) * (1:Nx);
  ys = (1:Ny)' * ones(1,Nx);
  xs = xs - halfX;
  ys = ys - halfY;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(Nx/2,Ny/2);
  radiusImg = img .* radiusMask;

  if mod( Nx, 2 )==0
    locs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * delta;
  else
    locs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * delta;
  end

  sinogram = zeros( nTheta, nDetectors );
  %parfor th=1:numel(thetas)
for th=1:numel(thetas)
    theta = thetas_deg(th);
    thisTrans_m = translations(th,:);
    thisTrans_pix = thisTrans_m / delta;
    translated = translateImg( radiusImg, thisTrans_pix );
    rotImg = imrotate( translated, theta, 'bilinear','crop' );
    sumResult = sum( rotImg, 1 ) * delta;

    interped = interp1( locs, sumResult, dLocs,'linear',0 );
    
    sinogram(th,:) = interped;
    if mod(th,10)==0
      disp(['ctRadon Theta: ', num2str(th), ' of ', ...
        num2str(numel(thetas)) ]);
    end;
  end

end

