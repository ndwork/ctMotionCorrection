
function [recon,costs] = ctCorrectForTranslation_LADMM( sinogram, ...
  nDetectors, detSize, thetas, translations_m, nCols, nRows, pixSize, ...
  varargin )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  defaultLambda = [];
  defaultMu = [];
  defaultR = [];
  p = inputParser;
  p.addOptional( 'lambda', defaultLambda, @isnumeric );
  p.addOptional( 'mu', defaultMu, @isnumeric );
  p.addOptional( 'radonMatrix', defaultR );
  p.parse( varargin{:} );
  lambda = p.Results.lambda;
  mu = p.Results.mu;
  R = p.Results.radonMatrix;
  
  gamma = 1d-6;   % Regularization parameter

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  if numel(R)==0
    R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
      detSize, thetas);
  end
  RT = transpose(R);

  translations_pix = translations_m / pixSize;

  applyE = @(u) RWithTranslation( u, translations_pix, nDetectors, RT );
  %applyE = @(u) radonWithTranslation( u, pixSize, nDetectors, ...
  %  detSize, thetas, translations_m );

  applyET = @(u) RTWithTranslation( u, translations_pix, nCols, RT );
  %cx = 0;  cy = 0;
  %applyET = @(u) backprojectionWithTranslation( u, thetas, ...
  %  detSize, cx, cy, nCols, nRows, pixSize, translations_m );
  %applyET = @(u) radonWithTranslationAdjoint( u, thetas, ...
  %  detSize, cx, cy, nCols, nRows, pixSize, translations_m );

  maxIters = 1000;
  x0 = rand( nRows, nCols );
  [nrmK, ~] = estimateNormKByPowerIteration( ...
    applyE, applyET, applyD1, applyD1T, applyD2, applyD2T, maxIters, x0 );

  if numel( lambda ) == 0 && numel( mu ) == 0
    minLambda = 1e-5; 
    maxLambda = 1e5;
    [lambda, mu] = findGoodStepSizes_LADMM( minLambda, maxLambda, nrmK, ...
      sinogram, nDetectors, detSize, thetas, translations_m, nCols, ...
      nRows, pixSize );
  elseif numel( lambda ) == 0
    mu = lambda/(nrmK*nrmK) * 0.9999;
  end
  %if u > lambda / (nrmK*nrmK)
  %  error('Improperly chosen step sizes');
  %end
  
  nThetas = numel( thetas );
  x = zeros( nRows, nCols );
  xBarE = zeros( nThetas, nDetectors );
  xBarD1 = zeros( nRows, nCols );
  xBarD2 = zeros( nRows, nCols );
  zE = zeros( nThetas, nDetectors );
  zD1 = zeros( nRows, nCols );
  zD2 = zeros( nRows, nCols );

  nIter = 1000;
  costs = zeros(nIter,1);
  minCost = 99999;
%reconH = figure;
  for i=1:nIter
    if mod(i,50)==0
      disp(['Working on iteration ', num2str(i), ' of ', num2str(nIter)]);
      %figure(reconH);  imshow( imresize(x,10,'nearest'), [] );
      %title(['Iteration ', num2str(i)]);  drawnow;
    end

    % Update x
    Ex = applyE(x);
    D1x = applyD1(x);
    D2x = applyD2(x);
    ET1 = applyET( Ex - zE + xBarE );
    ET2 = applyD1T( D1x - zD1 + xBarD1 );
    ET3 = applyD2T( D2x - zD2 + xBarD2 );
    argX = x - mu/lambda * ( ET1 + ET2 + ET3 );
    x = max( argX, 0 );

    % Update z
    Ex = applyE(x);
    D1x = applyD1(x);
    D2x = applyD2(x);
    argZE = Ex + xBarE;
    argZD1 = D1x + xBarD1;
    argZD2 = D2x + xBarD2;
    zE = ( lambda*sinogram + argZE ) / (lambda + 1);
    zD1 = softThresh( argZD1, lambda*gamma );
    zD2 = softThresh( argZD2, lambda*gamma );

    % Update xBar
    xBarE = xBarE + Ex - zE;
    xBarD1 = xBarD1 + D1x - zD1;
    xBarD2 = xBarD2 + D2x - zD2;

    % Store cost
    costs(i) = 0.5*norm( Ex(:) - sinogram(:), 2 )^2 + ...
      gamma * norm( D1x(:), 1 ) + gamma * norm( D2x(:), 1 );
    if costs(i) < minCost
      bestX = x;
      minCost = costs(i);
    end
  end
%close( reconH );

  recon = bestX;
end
