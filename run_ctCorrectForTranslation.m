
function run_ctCorrectForTranslation
  clear; close all;
  addpath(genpath('.'));

  % Reconstruction parameters
  method = 'PC';    % Options: GD, PC, LADMM
  cy = 0;   nRows=32;
  cx = 0;   nCols=32;
  pixSize = 0.001; % meters / pixel

  im = phantom();
  im = imresize( im, [nCols nRows], 'bilinear' );

  figure; imshow( imresize(im,10,'nearest'), [] );  title('original');

  detSize = 0.001;
  dTheta = 2 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);

  nDetectors = nCols*2;

  nonzeroTranslations = 1;
  if nonzeroTranslations > 0
    maxVerticalShift = 0.01; % in meters
    maxHorizontalShift = 0.02; % in meters
    translations = zeros( nThetas, 2 );
    translations(:,1) = linspace(0,maxVerticalShift,nThetas);
    translations(:,2) = linspace(0,maxHorizontalShift,nThetas);
  else
    maxVerticalShift = 0; % in meters
    maxHorizontalShift = 0; % in meters
    translations = zeros( nThetas, 2 );
  end

  im = padImgForRadon( im, maxHorizontalShift, maxVerticalShift, ...
      pixSize );
  [nRows,nCols] = size(im);

  sinogram = radonWithTranslation( im, pixSize, nDetectors, detSize, ...
   thetas, translations );

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

