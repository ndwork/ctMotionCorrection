
function run_ctCorrectForTranslation
  clear; close all;
  addpath(genpath('.'));

  % Reconstruction parameters
  method = 'LADMM';    % Options: GD, PC, LADMM
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
  translations(:,1) = linspace(0,maxVerticalShift,nThetas);
  translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

  % pad the phantom image so that it always stays in the field of view
  xShiftPix = maxHorizontalShift / pixSize;
  xPadding = zeros(size(im,1),ceil(abs(xShiftPix)));
  im = [xPadding im xPadding];
  yShiftPix = maxVerticalShift / pixSize;
  yPadding = zeros(ceil(abs(yShiftPix)),size(im,2));
  im = [ yPadding; im; yPadding ];

  sinogram = radonWithTranslation( im, pixSize, nDetectors, detSize, ...
   thetas, translations );
%  save('phSinogram.mat','sinogram')
%   load 'phSinogram.mat'

  %profile on
  tic;
  [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, detSize, ...
    thetas, translations, nCols, nRows, pixSize, 'method', method );
  timeTaken = toc;
  %profile off
  %profile viewer

  disp(['Time taken: ', num2str(timeTaken)]);
  figure; imshow( imresize(recon,10,'nearest'), [] );
  title('Reconstructed image');
  
  figure; plot( costs, 'LineWidth', 2 );
  xlabel('Iteration'); ylabel('Cost Function');

end

