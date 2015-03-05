function out = radonWithTranslationAdjoint( sinogram, thetas, detSize, ...
  cx, cy, Nx, Ny, pixSize, translations_m )
  [nThetas, nDetectors] = size(sinogram);

  dOffset=0;  % center channel offset

  nThetas = numel(thetas);
  dLocs = ( [0:nDetectors-1] - 0.5*(nDetectors-1) ) * detSize - dOffset;
  thetas_deg = thetas * 180/pi;

  if mod( Nx, 2 )==0
    pixLocsX = ( (0:Nx-1) - 0.5*Nx + 0.5 ) * pixSize;
  else
    pixLocsX = ( (0:Nx-1) - floor(0.5*Nx) ) * pixSize;
  end

  translations_pix = translations_m / pixSize;

  xs = ones(Ny,1) * (1:Nx);
  ys = (1:Ny)' * ones(1,Nx);
  halfY = Ny/2;  ys = ys - halfY;
  halfX = Nx/2;  xs = xs - halfX;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(Nx/2,Ny/2);

  M = makeLinearInterpMatrix( pixLocsX, dLocs );

  out = zeros(Ny,Nx);
  %parfor th = 1:nThetas
for th = 1:nThetas
    %extrapVal = 0;
    %interped = interp1( dLocs, sinogram(th,:), pixLocsX, 'linear', ...
    %  extrapVal );
    interped = M' * sinogram(th,:)';

    smear = ones(Ny,1) * interped';

    theta = thetas_deg(th);
    rotImg = imrotate( smear, -theta, 'bilinear', 'crop' );

    thisTrans_pix = translations_pix(th,:);
    translated = translateImg( rotImg, -thisTrans_pix );

    masked = translated .* radiusMask;
    
    out = out + masked;
  end

end
