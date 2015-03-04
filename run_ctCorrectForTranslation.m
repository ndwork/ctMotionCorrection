
function run_ctCorrectForTranslation
  clear; close all;
  addpath(genpath('.'));

  % Reconstruction parameters
  method = 'GD';    % Options: GD, PC, LADMM
  cx = 0;   nCols=32;
  cy = 0;   nRows=32;
  pixSize = 0.001; % meters / pixel

  im = phantom();
  im = imresize( im, [nCols nRows], 'bilinear' );

figure; imshow( imresize(im,10,'nearest'), [] );  title('original');

  detSize = 0.001;
  dTheta = 2 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);

  nDetectors = nCols*2;

  maxVerticalShift = 0.01; % in meters
  maxHorizontalShift = 0.02; % in meters
  translations = zeros( nThetas, 2 );
%   translations(:,1) = linspace(0,maxVerticalShift,nThetas);
%   translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

  % pad the phantom image so that it always stays in the field of view
  xShiftPix = maxHorizontalShift / pixSize;
  xPadding = zeros(size(im,1),ceil(abs(xShiftPix)));
  im = [xPadding im xPadding];
  yShiftPix = maxVerticalShift / pixSize;
  yPadding = zeros(ceil(abs(yShiftPix)),size(im,2));
  im = [ yPadding; im; yPadding ];

%   sinogram = radonWithTranslation( im, pixSize, nDetectors, detSize, ...
%    thetas, translations );
%   save('phSinogram.mat','sinogram')
  load 'phSinogram.mat'
  
%   maxIters = 1000;
%   x0 = rand( nRows, nCols );
%   [nrmK, lambdaVals] = estimateNormKByPowerIteration( ...
%    applyE, applyET, applyD1, applyD1T, applyD2, applyD2T, maxIters, x0 );
%   figure;  plot(lambdaVals);  title('Lambda v Iteration');
%   save( 'nrmK.mat', 'nrmK', 'lambdaVals' );
  load 'nrmK.mat';

%   minStep = 1e-5; 
%   maxStep = 1e5;
%   [optimalSigma, optimalTau] = findBestStepSizes(minStep,... 
%   maxStep, minStep, maxStep, nrmK, sinogram,...
%   nDetectors, detSize, thetas, translations, nCols, nRows, pixSize);
% %   load 'optimalSteps.mat'
%   save( 'optimalSteps.mat','optimalSigma', 'optimalTau' );


profile on
  tic;
  [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, detSize, ...
    thetas, translations, nCols, nRows, pixSize, nrmK, 'method', method );
%   [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, detSize, ...
%     thetas, translations, nCols, nRows, pixSize, nrmK, ...
%     optimalSigma, optimalTau );
  timeTaken = toc;
profile off
profile viewer

  disp(['Time taken: ', num2str(timeTaken)]);
  figure; imshow( recon, [] );  title('Reconstructed image');
  
  figure; plot( costs, 'LineWidth', 2 );
  xlabel('Iteration'); ylabel('Cost Function');
  
end



