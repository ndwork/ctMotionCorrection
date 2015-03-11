function [finaltranslation] = run_autofocus
  clear; close all;
  addpath(genpath('.'));

  % Reconstruction parameters
  method = 'PC';    % Options: GD, PC, LADMM
  cy = 0;   nRows=32;
  cx = 0;   nCols=32;
  pixSize = 0.001; % meters / pixel

  im = phantom();
  im = imresize( im, [nCols nRows], 'bilinear' );

  figure; imshow( imresize(im,10,'nearest'), [] );
  title('original'); drawnow;

  detSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);

  nDetectors = nCols*2;

  nonzeroTranslations = 1;
  if nonzeroTranslations > 0
    maxVerticalShift = 0.02; % in meters
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

  minX = 0;
  minY = 0;
  maxX = 0.02;
  maxY = 0.02;
  dx = 0.02;
  dy = 0.02;

  [finaltranslation] = autofocusWithTranslation(sinogram, minX, maxX, dx, minY, maxY, dy, nRows, nCols, nDetectors, pixSize, detSize, thetas, method);

end