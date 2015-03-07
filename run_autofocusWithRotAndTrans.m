function [finaltranslation] = run_autofocusWithRotAndTrans
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
  maxHorizontalShift = -0.1;
  maxVerticalShift = 0.1;
  
  translations = zeros( nThetas, 2 );
  translations(:,1) = linspace(0,maxVerticalShift,nThetas);
  translations(:,2) = linspace(0,maxHorizontalShift,nThetas);
  
  %Rotation of Image
  maxRotation = 0 * pi/180;
  rotations = linspace(0,maxRotation,nThetas);
  
  im = phantom();
  im = imresize( im, [nCols nRows], 'bilinear' );

  figure; imshow( imresize(im,10,'nearest'), [] );  title('original');
  
  im = padImgForRadon( im, maxHorizontalShift, maxVerticalShift, ...
      pixSize );
  
  sinogram = radonWithRotAndTrans( im, pixSize, nDetectors, detSize, ...
   thetas, rotations, translations );

  minX = -0.1;
  minY = 0;
  maxX = 0.1;
  maxY = 0;
  maxRot = 100*pi/180;
  dx = .1;
  dy = .1;
  dr = 50*pi/180;

  [finaltranslation] = autofocusWithRotAndTrans(sinogram, minX, maxX, dx, minY, maxY, dy, maxRot, dr, nRows, nCols, nDetectors, pixSize, detSize, thetas, method);
    
end