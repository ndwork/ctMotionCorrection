function bp = backprojectionWithTranslation( sinogram, thetas, detSize, ...
  cx, cy, Nx, Ny, pixSize, translations )

  [nThetas, nDetectors] = size(sinogram);

  dOffset = 0;   % detector center offset
  dLocs = ( (0:nDetectors-1) - 0.5*(nDetectors-1) ) * detSize - dOffset;

  % Make arrays of x and y positions of each pixel
  if mod( Nx, 2 )==0
    lineXs = ( (0:Nx-1) - 0.5*Nx + 0.5 ) * pixSize + cx;
  else
    lineXs = ( (0:Nx-1) - floor(0.5*Nx) ) * pixSize + cx;
  end
  if mod( Ny, 2 )==0
    lineYs = ( (0:Ny-1) - 0.5*Ny + 0.5 ) * pixSize + cy;
  else
    lineYs = ( (0:Ny-1) - floor(0.5*Ny) ) * pixSize + cy;
  end
  xs = ones(Ny,1) * lineXs;
  ys = lineYs' * ones(1,Nx);
  xs=xs(:) + cx;
  ys=ys(:) + cy;

  angles = atan2(ys,xs);
  pixDs = sqrt( xs.*xs + ys.*ys );

  xs = ones(Ny,1) * (1:Nx);
  ys = (1:Ny)' * ones(1,Nx);
  halfY = Ny/2;  ys = ys - halfY;
  halfX = Nx/2;  xs = xs - halfX;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(Nx/2,Ny/2);

  bp = zeros(Ny,Nx);
  parfor thIndx = 1:nThetas
    theta = thetas( thIndx );

    projections = pixDs .* cos( angles - theta );
    interped = interp1( dLocs, sinogram(thIndx,:), projections, ...
      'linear', 'extrap') * pixSize;
    interpedImg = reshape( interped, Ny, Nx );
    
    masked = interpedImg .* radiusMask;

    thisTrans_m = translations(thIndx,:);
    thisTrans_pix = thisTrans_m ./ [pixSize,pixSize];
    translated = translateImg( masked, -thisTrans_pix );

    bp = bp + translated;
  end

end
