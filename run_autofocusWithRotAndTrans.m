function [finaltranslation] = run_autofocusWithTranslation()
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

  nonzeroRotations = 1;
  if nonzeroRotations > 0
    maxRotation = 6 * pi/180;
    rotations = linspace(0,maxRotation,nThetas);
  else
    maxRotation = 0;
    rotations = zeros(nThetas,1);
  end
  
  im = padImgForRadon( im, maxHorizontalShift*1.5, ...
    maxVerticalShift*1.5, pixSize );
  [nRows,nCols] = size(im);

  sinogram = radonWithRotAndTrans( im, pixSize, nDetectors, detSize, ...
   thetas, rotations, translations );

  maxX = maxHorizontalShift*1.5;
  maxY = maxVerticalShift*1.5;
  dx = maxHorizontalShift/2;
  dy = maxVerticalShift/2;

  [finalTranslation,metric] = autofocusWithRotAndTrans(sinogram, ...
    maxX, dx, maxY, dy, maxRot, dr, nRows, nCols, pixSize, detSize, ...
    thetas, method);

  translations = zeros( nThetas, 2 );
  translations(:,1) = linspace(0,finalTranslation(1),nThetas);
  translations(:,2) = linspace(0,finalTranslation(2),nThetas);
  autofocusRecon = ctCorrectForTranslation( sinogram, nDetectors, ...
    detSize, thetas, translations, nCols, nRows, pixSize, ...
    'method', method, 'radonMatrix', R );
  
  figure; imshow( imresize(autofocusRecon,10,'nearest'), [] );
  title('Autofocus Reconstruction');
  
  figure; imshow( imresize(metric,10,'nearest'), [] );
  title('Autofocus Metric');
end
