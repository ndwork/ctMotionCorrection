function [] = makeFigures()

  close all; clear

  %% cartoon explaining what modified radon transform is

  % this is in run_makeImagesForReport.m
  
  %% PC results - no translation, with translation, with translation and rotation

  % Reconstruction parameters
  cy = 0;   nRows=5;
  cx = 0;   nCols=5;
  pixSize = 0.001; % meters / pixel

  im = phantom();
  im = imresize( im, [nCols nRows], 'bilinear' );

  figure; imshow( imresize(im,10,'nearest'), [] );  title('original');

  detSize = 0.001;
  dTheta = 2 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);

  nDetectors = nCols*2;
  
  maxShift = [0     0; 
              0.01  0.02; 
              0.01  0.02]; % [maxVertical maxHorizontal]
  rotation = [0               0;
              0               0;
              -10 * pi/180    10 * pi/180]; % [maxRot minRot]
  imageNames = {'NoTransNoRot','TransAndNoRot','TransAndRot'};
  
  method = {'PC','LADMM'};    % Options: GD, PC, LADMM
  
  for m = 1:numel(method)
                
    for k = 1:numel(imageNames)

      translations = zeros( nThetas, 2 );
      translations(:,1) = linspace(0,maxShift(k,1),nThetas);
      translations(:,2) = linspace(0,maxShift(k,2),nThetas);

      rotations = linspace(rotation(k,1),rotation(k,2),nThetas);

      im = padImgForRadon( im, maxShift(k,2), maxShift(k,1), ...
        pixSize );
      [nRows,nCols] = size(im);

      sinogram = radonWithRotAndTrans( im, pixSize, nDetectors, detSize, ...
       thetas, rotations, translations );

      if (rotation(k,1) == 0) && (rotation(k,2) == 0)
        [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, ...
          detSize, thetas, translations, nCols, nRows, pixSize, ...
          'method', char(method(m)) );
      else
        [recon,costs] = ctCorrectForRotAndTrans( sinogram, nDetectors, ...
          detSize, thetas, rotations, translations, nCols, nRows, pixSize, ...
          'method', char(method(m)) );
        
      end
      
      figure; imshow( imresize(recon,10,'nearest'), [] );
      fileName = [char(pwd) '/figures/' 'img' char(method(m)) ...
        char(imageNames(k)) '.eps'];
      saveas(gcf,fileName,'epsc');
      title('Reconstructed image');

      figure; set(gca,'FontSize',14) 
      plot( costs, 'LineWidth', 2 );
      xlabel('Iteration','FontSize',14); 
      ylabel('Cost Function','FontSize',14);
      fileName = [char(pwd) '/figures/' 'costs' char(method(m)) ...
        char(imageNames(k)) '.eps'];
      saveas(gcf,fileName,'epsc');

    end
  end
              

  
  
  %% LADMM results - no translation, with translation, with translation and rotation

  
  %% GD results - no translation, with translation, with translation and rotation

  
  
  %% autofocus result - no translation, with translation, with rotation and translation

  % this is in makeBlurryPhantoms.m
  images = makeBlurryPhantoms_copy();
  image1 = images(:,:,1);
  figure; imshow(image1,[])
  image2 = images(:,:,2);
  figure; imshow(image2,[])
  image3 = images(:,:,3);
  figure; imshow(image3,[])
  
  %% compare to inverse radon transform - no translation, with translation, with rotation and translation



end