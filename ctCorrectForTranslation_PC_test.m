
function [recon,costs] = ctCorrectForTranslation_PC_test( sinogram,...
  applyE, applyET, nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize, gamma,...
  sigma, tau )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram
  
  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));
    
%     maxIters = 1000;
%   x0 = rand( nRows, nCols );
%   [nrmK, lambdaVals] = estimateNormKByPowerIteration( ...
%    applyE, applyET, applyD1, applyD1T, applyD2, applyD2T, maxIters, x0 );
%   figure;  plot(lambdaVals);  title('Lambda v Iteration');
%   save( 'nrmK_deblur.mat', 'nrmK', 'lambdaVals' );
  load 'nrmK_deblur.mat';
  
  if nargin < 12
%     sigma = 1/nrmK;
%     tau = 0.99*(1/nrmK);

      tau = 10;
      sigma = 1/(tau*nrmK^2);
  end;

  x = zeros( nRows, nCols );
  xBar = zeros( nRows, nCols );
  yE = zeros( size(sinogram,1), size(sinogram,2) );
  yD1 = zeros( nRows, nCols );
  yD2 = zeros( nRows, nCols );
% 
%   x = sinogram;
%   xBar = x;
%   yE = applyE(x);
%   yD1 = applyD1(x);
%   yD2 = applyD2(x);

  if sigma*tau > 1 / (nrmK*nrmK)
    error('Improperly chosen step sizes');
  end

  alpha = 1;
  nIter = 5000;
  costs = zeros(nIter,1);
  minCost = 9999;  bestX = x;
reconH = figure;
  for i=1:nIter
    if mod(i,5)==0
      disp(['Working on iteration ', num2str(i), ' of ', num2str(nIter)]);
%       figure(reconH);  imshow( imresize(x,10,'nearest'), [] );
      figure(reconH);  imshow(x, [] );
      title(['Iteration ', num2str(i)]);  drawnow;
    end

    % Update y
    ExBar = applyE( xBar );     tmpE = yE + sigma * ExBar;
    D1xBar = applyD1( xBar );   tmpD1 = yD1 + sigma * D1xBar;
    D2xBar = applyD2( xBar );   tmpD2 = yD2 + sigma * D2xBar;

    yE = ( tmpE - sigma*sinogram ) / ( sigma + 1 );
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

