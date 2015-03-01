
function [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize, nrmK, ...
  sigma, tau )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  if nargin < 10, sigma=1; end;

  %gamma = 1d-5;   % Regularization parameter
  gamma = 0;

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  applyE = @(u) radonWithTranslation( u, pixSize, nDetectors, ...
    detSize, thetas, translations );

  cx = 0;  cy = 0;
  applyET = @(u) backprojectionWithTranslation( u, thetas, ...
    detSize, cx, cy, nCols, nRows, pixSize, translations );

  nThetas = numel( thetas );
  x = zeros( nRows, nCols );
  xBar = zeros( nRows, nCols );
  yE = zeros( nThetas, nDetectors );
  yD1 = zeros( nRows, nCols );
  yD2 = zeros( nRows, nCols );

  if nargin < 10, tau = 1 / (nrmK*nrmK*sigma) * 0.999; end;

  if sigma*tau > 1 / (nrmK*nrmK)
    error('Improperly chosen step sizes');
  end

  alpha = 1;
  nIter = 1000;
  costs = zeros(nIter,1);
  minCost = 9999;  bestX = x;
reconH = figure;
  for i=1:nIter
    if mod(i,2)==0
      disp(['Working on iteration ', num2str(i), ' of ', num2str(nIter)]);
      figure(reconH);  imshow( imresize(x,10,'nearest'), [] );  drawnow;
    end

    % Update y
    ExBar = applyE( xBar );     tmpE = yE + sigma * ExBar;
    D1xBar = applyD1( xBar );   tmpD1 = yD1 + sigma * D1xBar;
    D2xBar = applyD2( xBar );   tmpD2 = yD2 + sigma * D2xBar;

 
    %tmpBefore = 0.5*yE(:)'*yE(:) + yE(:)'*sinogram(:) + 1/(2*sigma)*norm(yE(:)-tmpE(:),2)^2;
    yE = ( tmpE - sigma*sinogram ) / ( sigma + 1 );
    %tmpAfter  = 0.5*yE(:)'*yE(:) + yE(:)'*sinogram(:) + 1/(2*sigma)*norm(yE(:)-tmpE(:),2)^2;
    %disp(['Should be negative: ', num2str(tmpAfter-tmpBefore)]);
    yD1 = min( tmpD1, gamma );
    yD1 = max( yD1, -gamma );
    yD2 = min( tmpD2, gamma );
    yD2 = max( yD2, -gamma );
    
    % Store cost
    costs(i) = 0.5*norm( ExBar(:) - sinogram(:), 2 )^2 + ...
      gamma * norm( D1xBar(:), 1 ) + gamma * norm( D2xBar(:), 1 );
    if costs(i) < minCost
      minCost = costs(i);
      bestX = x;
    end

    % Update x
    lastX = x;
    ETyE = applyET(yE);
    D1TyD1 = applyD1T(yD1);
    D2TyD2 = applyD2T(yD2);
    tmp = x - tau * ( ETyE + D1TyD1 + D2TyD2 );
    x = max( tmp, 0 );
    
    % Update xBar
    xBar = x + alpha * ( x - lastX );
  end

  recon = bestX;
end

