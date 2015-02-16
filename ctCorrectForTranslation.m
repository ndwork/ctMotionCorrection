
function recon = ctCorrectForTranslation( sinogram, delta, nDetectors, ...
  dSize, thetas, translations )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  sizeSino = size(sinogram);
  nRows = sizeSino(2);  %num rows of reconstructed image
  nCols = sizeSino(2);  %num cols of reconstructed image

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1,:), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  applyE = @(u) radonWithTranslation( u, delta, nDetectors, ...
    dSize, thetas, translations );
  
  cx = 0;  cy = 0;
  applyET = @(u) applyETrans( u, thetas, dSize, cx, cy, nCols, nRows, ...
    dSize, dSize, translations );


  maxIters = 1000;
  x0 = rand( sizeSino(2) );
  [nrmK, lambdaVals] = estimateNormByPowerIterationED( ...
    applyE, applyET, applyD1, applyD1T, applyD2, applyD2T, maxIters, x0 );
  figure;  plot(lambdaVals);  title('Lambda v Iteration');

save( 'nrmK.mat', 'nrmK', 'lambdaVals' );
  
  xE = zeros( nRows, nCols );
  xD1 = zeros( nRows, nCols );
  xD2 = zeros( nRows, nCols );
  xBarE = zeros( sizeSino );
  xBarD1 = zeros( sizeSino );
  xBarD2 = zeros( sizeSino );
  yE = zeros( sizeSino );
  yD1 = zeros( nRows, nCols );
  yD2 = zeros( nRows, nCols );


  nIter = 1000;

  if sigma*tau > 1 / (nrmK*nrmK)
    error('Improperly chosen step sizes');
  end

  for i=1:nIter

  end



end


function bp = applyETrans( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy, translations )

  [nThetas, nDetectors] = size(sinogram);

  dOffset = 0;   % detector center offset
  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  % Make arrays of x and y positions of each pixel
  if mod( Nx, 2 )==0
    lineXs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * dx + cx;
  else
    lineXs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * dx + cx;
  end
  if mod( Ny, 2 )==0
    lineYs = ( ([0:Ny-1]) - 0.5*Ny + 0.5 ) * dy + cy;
  else
    lineYs = ( ([0:Ny-1]) - floor(0.5*Ny) ) * dy + cy;
  end
  xs = ones(Ny,1) * lineXs;
  ys = lineYs' * ones(1,Nx);
  xs=xs(:) + cx;
  ys=ys(:) + cy;

  angles = atan2(ys,xs);
  pixDs = sqrt( xs.*xs + ys.*ys );

  bp = zeros(Ny,Nx);
  parfor thIndx = 1:nThetas
    theta = thetas( thIndx );

    projections = pixDs .* cos( angles - theta );
    interped = interp1( dLocs, sinogram(thIndx,:), projections, ...
      'linear', 0);

    interpedImg = reshape( interped, Ny, Nx );

    thisTrans_m = translations(thIndx,:);
    thisTrans_pix = thisTrans_m ./ [dy,dx];
    translated = translateImg( interpedImg, -thisTrans_pix );

    bp = bp + translated;
    if mod(thIndx,10)==0
      disp([ 'applyETrans Theta Indx / Theta: ', ...
        num2str(thIndx), '/', num2str(theta) ]);
    end
  end
end



