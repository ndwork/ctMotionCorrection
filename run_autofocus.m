function [finaltranslation] = run_autofocus
  method = 'PC';    % Options: GD, PC, LADMM
  nRows=32;
  nCols=32;
  nDetectors = nCols*2;
  pixSize = 0.001; % meters / pixel
  detSize = 0.001;
  
  dTheta = 2 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);
  
  %Translation of Image
  maxHorizontalShift = -0.01;
  maxVerticalShift = 0.01;
  
  translations = zeros( nThetas, 2 );
  translations(:,1) = linspace(0,maxVerticalShift,nThetas);
  translations(:,2) = linspace(0,maxHorizontalShift,nThetas);
  
  im = phantom();
  im = imresize( im, [nCols nRows], 'bilinear' );

  figure; imshow( imresize(im,10,'nearest'), [] );  title('original');
  
  im = padImgForRadon( im, maxHorizontalShift, maxVerticalShift, ...
      pixSize );
  
  sinogram = radonWithTranslation( im, pixSize, nDetectors, detSize, ...
   thetas, translations );

  minX = -0.01;
  minY = -0.01;
  maxX = 0.01;
  maxY = 0.01;
  dx = .01;
  dy = .01;

  [finaltranslation] = autofocus(sinogram, minX, maxX, dx, minY, maxY, dy, nRows, nCols, nDetectors, pixSize, detSize,thetas, method);

end