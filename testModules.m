function [] = testModules()
  close all; clear;
  profile on;


  %% Test making linear interpolation matrix
  x = rand(20,1);
  xLocs = [1:1:20]';
  interpLocs = rand(20,1)*20;

  interpMatlab = interp1(xLocs,x,interpLocs,'linear',0);
  interpMatrix = makeLinearInterpMatrix( xLocs, interpLocs );
  interpUs = interpMatrix*x;

  error = max(abs(interpMatlab - interpUs));
  if error < 1d-12
    disp('Test of linear interpolation matrix: Passed');
  else
    disp(['Test of linear interpolation matrix: Failed with error ',...
      num2str(error)]);
  end

  %% Test image translation function
  x = phantom();
  trans = [-12 -3];
  trueTrans = circshift( x, trans );
  xTrans = translateImg(x,trans);
  error = max( abs( trueTrans(:) - xTrans(:) ) );
  if error < 1d-12
    disp('Test of translate image: Passed');
  else
    disp(['Test of translate image: Failed with error ', num2str(error)]);
  end

  %% Test adjoint of D1
  nDetectors = 500;
  nRows = nDetectors;
  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1)); 
  [outD1, adjointError] = testAdjoint(applyD1,applyD1T,nDetectors,nDetectors,nDetectors,nDetectors);
  if outD1 == 1;
    disp('test adjoint of D1: Passed')
  else
    disp(['test adjoint of D1: Failed with adjoint error: ', ...
      num2str(adjointError)]')
  end

  %% Test adjoint of D2
  nDetectors = 500;
  nCols = nDetectors;
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));
  [outD2, adjointError] = testAdjoint(applyD2,applyD2T,nDetectors,nDetectors,nDetectors,nDetectors);
  if outD2 == 1;
    disp('test adjoint of D2: Passed')
  else
    disp(['test adjoint of D2: Failed with adjoint error: ', ...
      num2str(adjointError)]')
  end

  %% test adjoint of Translation
  outT = testAdjoint_translateImg(100,150,100,150);
  if outT == 1;
    disp('test adjoint of translateImg: Passed')
  else
    disp('test adjoint of translateImg: Failed')
  end

  %% simple test adjoint of R (E without translation)
  % this does not take into account the fact that the image of the Radon
  % transforms is sinograms
  outR = simpleTestAdjointR();

  if outR == 1;
    disp('simple test adjoint of R: Passed')
  else
    disp('simple test adjoint of R: Failed')
  end

  %% test adjoint of R (E without translation)
  nDetectors = 500;
  detectorSize = 1;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);
  translations = zeros( nThetas, 2 );
  sizeSino = [nThetas nDetectors];
  nRows = sizeSino(2)/2;  %num rows of reconstructed image
  nCols = sizeSino(2)/2;  %num cols of reconstructed image
  pixelSize = 1;

  applyR = @(u) radonWithTranslation( u, pixelSize, nDetectors, ...
    detectorSize, thetas, translations );

  cx = 0;  cy = 0;
  applyRT = @(u) radonWithTranslationAdjoint( u, thetas, detectorSize, ...
    cx, cy, nRows, nCols, pixelSize, translations );

  [outR,adjointError] = testAdjointRadon(applyR,applyRT,nRows,nCols);

  if outR == 1;
    disp('Test adjoint of R: Passed')
  else
    disp(['Test adjoint of R: Failed with adjoint error: ', ...
      num2str(adjointError)])
  end

  %% test adjoint of E
  nDetectors = 500;
  detectorSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);
  maxVerticalShift = 0.012;
  maxHorizontalShift = -0.028;
  translations = zeros( nThetas, 2 );
  translations(:,1) = linspace(0,maxVerticalShift,nThetas);
  translations(:,2) = linspace(0,maxHorizontalShift,nThetas);
  sizeSino = [nThetas nDetectors];
  nRows = sizeSino(2)/2;  %num rows of reconstructed image
  nCols = sizeSino(2)/2;  %num cols of reconstructed image
  pixelSize = 0.001;

  applyE = @(u) radonWithTranslation( u, pixelSize, nDetectors, ...
    detectorSize, thetas, translations );

  cx = 0;  cy = 0;
  %applyET = @(u) backprojectionWithTranslation( u, thetas, detectorSize, ...
  %  cx, cy, nRows, nCols, pixelSize, translations );
  applyET = @(u) radonWithTranslationAdjoint( u, thetas, detectorSize, ...
    cx, cy, nRows, nCols, pixelSize, translations );

  [outR,adjointError] = testAdjointRadon(applyE,applyET,nRows,nCols);

  if outR == 1;
    disp('Test adjoint of E: Passed')
  else
    disp(['Test adjoint of E: Failed with adjoint error: ', ...
      num2str(adjointError)])
  end

  %% Test estimating the norm of K by power iteration
  sizeX = 20;
  sizeY = 20;
  
  test = rand(sizeX,sizeY);
  applyA = @(u) test*u;
  applyATrans = @(u) test'*u;
  
  power = estimateNormByPowerIteration(applyA,applyATrans,...
    rand(sizeX,sizeY));
  power2 = normest(test);
  
  error = abs(power - power2);
  if error < 1e-9
    disp('Test of norm estimation: Passed')
  else
    disp(['Test of norm estimation: Failed with error: ' ...
      num2str(error)])
  end

    %% Test Radon Matrix
  im = phantom();
  im = imresize( im, [32 32] );
  [nCols, nRows] = size(im);
  pixSize = 1;
  nDetectors = 100;
  detSize = 1;
  dTheta = 2 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nthetas = numel(thetas);

  if exist('RadonMatrix.mat', 'file') == 2
    load 'RadonMatrix.mat';
  else
    R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
      detSize, thetas);
    save( 'RadonMatrix.mat', 'R');
  end

  sino1 = R * im(:);
  sino1 = reshape( sino1, [nthetas nDetectors] );

  sino2 = radonTransform( im, pixSize, nDetectors, ...
    detSize, thetas );

  error = max( abs( sino1(:) - sino2(:) ) );
  if error < 1d-12
    disp('Test of makeRadonMatrix: Passed');
  else
    disp(['Test of makeRadonMatrix: Failed with error ',...
      num2str(error)]);
  end

  %% Cleanup
  profile off;
  profile viewer;
end

function [out] = testAdjoint_translateImg(Mx,Nx,My,Ny)
  t1 = [1.5 2.2];
  t2 = -t1;

  apply = @(u) translateImg(u,t1);
  applyTrans = @(u) translateImg(u,t2);

  x = rand(Mx,Nx);
  y = rand(My,Ny);

  x( :, 1:ceil(t1(2)) ) = 0;
  x( :, end-ceil(t1(2))+1:end) = 0;
  x(1:ceil(t1(1)), :) = 0 ;
  x(end-ceil(t1(1))+1:end, :) = 0;

  y( :, 1:ceil(t1(2)) ) = 0;
  y( :, end-ceil(t1(2))+1:end) = 0;
  y(1:ceil(t1(1)), :) = 0 ;
  y(end-ceil(t1(1))+1:end, :) = 0;

  tmp1 = sum(sum(apply(x).*y));
  tmp2 = sum(sum(x.*applyTrans(y)));

  if abs(tmp1 - tmp2) < 1e-10
    out = 1;
  else
    out = 0;
  end
end

function [out, error] = testAdjointRadon( applyR, applyRT, Mx, Nx )
  img = phantom();
  img = imresize(img,[Mx Nx]);
  Rimg = applyR(img);

  imgToMakeSino = imrotate(phantom(),90);
  sino = applyR(imgToMakeSino);
  RTsino = applyRT(sino);

  prod1 = Rimg .* sino;     innerProd1 = sum( prod1(:) );
  prod2 = img .* RTsino;    innerProd2 = sum( prod2(:) );

  error = abs(innerProd1 - innerProd2) / min(innerProd1,innerProd2);

  if error < 1e-10
    out = 1;
  else
    out = 0;
  end

end

function [outR] = simpleTestAdjointR()

  nDetectors = 4;
  dSize = 0.1;
  thetas = [0 -pi/2];
  nThetas = numel(thetas);
  delta = 0.1;
  translations = zeros(nThetas,2);

  x = [1 2 3 4; 5 6 7 8; 9 10 1 2; 3 4 5 6];
  y = [10 9 8 7; 6 5 4 3];

  Rx = radonWithTranslation( x, delta, nDetectors, ...
    dSize, thetas, translations );
  innerProd_Rx_y = sum(sum(Rx.*y));

  cx = 0;
  cy = 0;
  [Ny, Nx] = size(x);
  RTx = backprojectionWithTranslation( y, thetas, dSize, cx, cy, Nx, Ny, ...
    delta, translations );  
  innerProd_x_RTy = sum(sum(x.*RTx));

  if abs(innerProd_Rx_y - innerProd_x_RTy) < 1e-10
    outR = 1;
  else
    outR = 0;
  end
end