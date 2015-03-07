
function [recon,costs] = ctCorrectForTranslation_LADMM_test( sinogram, ...
  applyE, applyET, nDetectors, detSize, thetas, translations, nCols, ...
  nRows, pixSize, gamma, lambda, mu )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  if nargin < 12, lambda=1; end;

%   maxIters = 1000;
%   x0 = rand( nRows, nCols );
%   [nrmK, lambdaVals] = estimateNormKByPowerIteration( ...
%    applyE, applyET, applyD1, applyD1T, applyD2, applyD2T, maxIters, x0 );
%   figure;  plot(lambdaVals);  title('Lambda v Iteration');
%   save( 'nrmK.mat', 'nrmK', 'lambdaVals' );
  load 'nrmK_rand.mat';
  
  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  x = zeros( nRows, nCols );
  xBarE = zeros( size(sinogram,1), size(sinogram,2) );
  xBarD1 = zeros( nRows, nCols );
  xBarD2 = zeros( nRows, nCols );
  zE = zeros( size(sinogram,1), size(sinogram,2) );
  zD1 = zeros( nRows, nCols );
  zD2 = zeros( nRows, nCols );

  if nargin < 13, mu = lambda / nrmK^2 * 0.99; end;

  nIter = 1000;
  costs = zeros(nIter,1);
  minCost = 9999;  bestX = x;
reconH = figure;
  for i=1:nIter
    if mod(i,5)==0
      disp(['Working on iteration ', num2str(i), ' of ', num2str(nIter)]);
%       figure(reconH);  imshow( imresize(x,10,'nearest'), [] );
      figure(reconH);  imshow( x, [] );
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
    
    if costs(i) < minCost
      minCost = costs(i);
      bestX = x;
    end
    
    % Update z
    Ex = applyE(x);
    D1x = applyD1(x);
    D2x = applyD2(x);
    argZE = Ex + xBarE;
    argZD1 = D1x + xBarD1;
    argZD2 = D2x + xBarD2;
    zE = ( lambda*sinogram + argZE ) / (lambda + 1);
    zD1 = softThresh( argZD1, lambda );
    zD2 = softThresh( argZD2, lambda );

    % Update xBar
    xBarE = xBarE + Ex - zE;
    xBarD1 = xBarD1 + D1x - zD1;
    xBarD2 = xBarD2 + D2x - zD2;
  end

  recon = bestX;

end
