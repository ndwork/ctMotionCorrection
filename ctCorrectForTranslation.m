
function [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, detSize, ...
  thetas, translations, nCols, nRows, pixSize )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  sigma = 1;
  gamma = 1d-3;   % Regularization parameter
  
  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  applyE = @(u) radonWithTranslation( u, pixSize, nDetectors, ...
    detSize, thetas, translations );

  cx = 0;  cy = 0;
  applyET = @(u) backprojectionWithTranslation( sinogram, thetas, detSize, ...
    cx, cy, nCols, nRows, pixSize, translations );


  maxIters = 1000;
  x0 = rand( nRows, nCols );
  %[nrmK, lambdaVals] = estimateNormKByPowerIteration( ...
  %  applyE, applyET, applyD1, applyD1T, applyD2, applyD2T, maxIters, x0 );
  %figure;  plot(lambdaVals);  title('Lambda v Iteration');
load 'nrmK.mat';
%save( 'nrmK.mat', 'nrmK', 'lambdaVals' );

  nThetas = numel( thetas );
  x = zeros( nRows, nCols );
  xBar = zeros( nRows, nCols );
  yE = zeros( nThetas, nDetectors );
  yD1 = zeros( nRows, nCols );
  yD2 = zeros( nRows, nCols );

  tau = 1 / (nrmK*nrmK*sigma) * 0.999;
  
  if sigma*tau > 1 / (nrmK*nrmK)
    error('Improperly chosen step sizes');
  end

  theta = 1;
  nIter = 1000;
  costs = zeros(nIter,1);
  for i=1:nIter
    % Update y
    ExBar = applyE( xBar );     tmpE = yE + sigma * ExBar;
    D1xBar = applyD1( xBar );   tmpD1 = yD1 + sigma * D1xBar;
    D2xBar = applyD2( xBar );   tmpD2 = yD2 + sigma * D2xBar;

    costs(i) = 0.5*norm( ExBar(:) - sinogram(:), 2 ) + ...
      gamma * norm( D1xBar(:), 1 ) + gamma * norm( D2xBar(:), 1 );
    
    yE = ( tmpE - sigma*sinogram ) / ( sigma + 1 );
    yD1 = min( tmpD1, gamma );
    yD1 = max( yD1, -gamma );
    yD2 = min( tmpD2, gamma );
    yD2 = max( yD2, -gamma );

    % Update x
    lastX = x;
    tmp = x - tau * ( applyET(yE) + applyD1T(yD1) + applyD2T(yD2) );
    x = max( tmp, 0 );

    % Update xBar
    xBar = x + theta * ( x - lastX );
  end

  recon = x;
end

