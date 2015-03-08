
function [recon,costs] = ctCorrectForTranslation_LADMM( sinogram, ...
  nDetectors, detSize, thetas, translations_m, nCols, nRows, pixSize, ...
  lambda, mu )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

%   maxIters = 1000;
%   x0 = rand( nRows, nCols );
%   [nrmK, lambdaVals] = estimateNormKByPowerIteration( ...
%    applyE, applyET, applyD1, applyD1T, applyD2, applyD2T, maxIters, x0 );
%   figure;  plot(lambdaVals);  title('Lambda v Iteration');
%   save( 'nrmK.mat', 'nrmK', 'lambdaVals' );
  load 'nrmK.mat';

  if nargin < 9, lambda=1d2; end;
  
  gamma = 1d-5;   % Regularization parameter

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  %R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
  %  detSize, thetas);
load 'RadonMatrix.mat';
  RT = transpose(R);

  translations_pix = translations_m / pixSize;
  applyE = @(u) RWithTranslation( u, translations_pix, nDetectors, R );
  %applyE = @(u) radonWithTranslation( u, pixSize, nDetectors, ...
  %  detSize, thetas, translations_m );

  %cx = 0;  cy = 0;
  %applyET = @(u) backprojectionWithTranslation( u, thetas, ...
  %  detSize, cx, cy, nCols, nRows, pixSize, translations_m );
  %applyET = @(u) radonWithTranslationAdjoint( u, thetas, ...
  %  detSize, cx, cy, nCols, nRows, pixSize, translations_m );
  applyET = @(u) RTWithTranslation( u, translations_pix, nCols, RT );

  nThetas = numel( thetas );
  x = zeros( nRows, nCols );
  xBarE = zeros( nThetas, nDetectors );
  xBarD1 = zeros( nRows, nCols );
  xBarD2 = zeros( nRows, nCols );
  zE = zeros( nThetas, nDetectors );
  zD1 = zeros( nRows, nCols );
  zD2 = zeros( nRows, nCols );

  if nargin < 10, mu = lambda / nrmK^2 * 0.9; end;

  nIter = 10000;
  costs = zeros(nIter,1);
reconH = figure;
  for i=1:nIter
    if mod(i,5)==0
      disp(['Working on iteration ', num2str(i), ' of ', num2str(nIter)]);
      figure(reconH);  imshow( imresize(x,10,'nearest'), [] );
      title(['Iteration ', num2str(i)]);  drawnow;
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

    % Store cost
    costs(i) = 0.5*norm( Ex(:) - sinogram(:), 2 )^2 + ...
      gamma * norm( D1x(:), 1 ) + gamma * norm( D2x(:), 1 );
    
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
  end
close( reconH );

  recon = bestX;
end
