function [] = makeFigures()

  close all; clear
  runAutoFocus = 0;
  runPCandLADMM = 0;
  runFilteredBackProjection = 1;
  runBlurryPhantom = 0;
  runBlurryPhantomMetrics = 0;
  runAnimation = 0;
  
  
  %% PC, LADMM, and GD results
  % no translation, with translation, with translation and rotation
  
    close all;

    % Reconstruction parameters
    noiseLevel = 0;
    cy = 0;   nRows=32;
    cx = 0;   nCols=32;
    pixSize = 0.001; % meters / pixel

    img = phantom();
    img = imresize( img, [nCols nRows], 'bilinear' );

    figure; imshow( imresize(img,10,'nearest'), [] );  
    fileName = [char(pwd) '/figures/originalImage.eps'];
    saveas(gcf,fileName,'epsc');
    title('original');

    detSize = 0.001;
    dTheta = 1 * pi/180;
    thetas = 0:dTheta:pi-dTheta;
    nThetas = numel(thetas);

    nDetectors = nCols*2;

    maxShift = [0     0; 
                0.01  0.02; 
                0.01  0.02]; % [maxVertical maxHorizontal]
% maxShift = [0.01  0.02; 
%                 0.01  0.02]; % [maxVertical maxHorizontal]
    rotation = [0               0;
                0               0;
                -10 * pi/180    10 * pi/180]; % [maxRot minRot]
% rotation = [0               0;
%                 0               0]; % [maxRot minRot]
    imageNames = {'NoTransNoRot','TransAndNoRot','TransAndRot'};
%     imageNames = {'TransAndNoRot','TransAndRot'};

    method = {'PC','LADMM'};    % Options: GD, PC, LADMM
%     method = {'PC'};    % Options: GD, PC, LADMM
    tic;

    
    for m = 1:numel(method)
%     FBPrun = 0;

      for k = 1:numel(imageNames)

        translations = zeros( nThetas, 2 );
        translations(:,1) = linspace(0,maxShift(k,1),nThetas);
        translations(:,2) = linspace(0,maxShift(k,2),nThetas);

        rotations = linspace(rotation(k,1),rotation(k,2),nThetas);

        im = padImgForRadon( img, maxShift(k,2), maxShift(k,1), ...
          pixSize );
        [nRows,nCols] = size(im);

        sinogram = radonWithRotAndTrans( im, pixSize, nDetectors, detSize, ...
           thetas, rotations, translations );

        noise = noiseLevel*randn(size(sinogram,1),size(sinogram,2));
        sinogram = sinogram + noise;
       
       if runFilteredBackProjection == 1 && m == 1
        recon = filteredBackProjection(sinogram, thetas, detSize, ...
            nCols, nRows, pixSize);
          
        figure; imshow( imresize(recon,10,'nearest'), [] );
        fileName = [char(pwd) '/figures/' 'imgFBP' ...
          char(imageNames(k)) '.eps'];
        saveas(gcf,fileName,'epsc');
        title('Reconstructed image');
         
%          FBPrun = 1;
       end
       
       if runPCandLADMM == 1;

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
    end
    timeTaken = toc;
    display(['Time taken: ', num2str(timeTaken)])
 
  %% cartoon explaining what modified radon transform is
  if runAnimation == 1;
  
    close all;

    img = phantom();
    translations = [0.5,0]; %meters
%     translations = [0,0]; %meters
    pixsize = 0.01;
    translations = translations/pixsize; %pixels
    rotations = -120; %degrees
%     rotations = 0; %degrees
    img = imrotate(img,-rotations/2,'crop');
    img = padImgForRadon(img,50,50,1);
  %     figure;
  %     imshow(img,[])
    numFrames = 4;
    images =  makeImagesForReport(img,rotations,translations,numFrames);

    for k = 1:numFrames
      thisImage = images(:,:,k);
      figure; imshow(thisImage, [])
      fileName = [char(pwd) '/figures/' 'phantomRollingImageTransAndRot' ...
        num2str(k) '.eps'];
      saveas(gcf,fileName,'epsc'); 
    end  
  end
  
  %% autofocus result - no translation, with translation, with rotation and translation
  if runAutoFocus == 1;
    images = 1;
    
  end
  
  %% make blurry phantom to demonstrate image metrics
  if runBlurryPhantom == 1;
  
    close all;

    horizontaltrans = [5 2 0];
    images = makeBlurryPhantoms(horizontaltrans);

    for k = 1:numel(horizontaltrans)
      thisImage = images(:,:,k);
      figure; imshow(thisImage, [])
      fileName = [char(pwd) '/figures/' 'blurryPhantomWithTrans' ...
        num2str(horizontaltrans(k)) '.eps'];
      saveas(gcf,fileName,'epsc');
    end
  end
  %% plot image metrics for blurry phantom
    % plot image metrics
    
  if runBlurryPhantomMetrics == 1
    
    [x, metric_normgrad, metric_laplacian, metric_histogram, ...
    metric_variance] = makeMetricImages();

    lw = 2;
    fs = 18; % font size for labels
    ts = 16; % font size for tick marks
    ys = 6; % number of y tick marks to show

    figure;
    set(gca,'FontSize',ts)
    plot(x,metric_normgrad,'LineWidth',lw);
    xlabel('Image Shift (pixels)','FontSize',fs)
    ylabel('Normalized Gradient Squared','FontSize',fs)
    set(gca,'XTick',linspace(-5,5,3))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),ys))
    set(gca,'XMinorTick','on','YMinorTick','on')
    fileName = [char(pwd) '/figures/' 'metric_normgrad.eps'];
    saveas(gcf,fileName,'epsc');

    figure;
    set(gca,'FontSize',ts)
    plot(x,metric_laplacian,'LineWidth',lw);
    xlabel('Image Shift (pixels)','FontSize',fs)
    ylabel('Laplacian','FontSize',fs)
    set(gca,'XTick',linspace(-5,5,3))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),ys))
    set(gca,'XMinorTick','on','YMinorTick','on')
    fileName = [char(pwd) '/figures/' 'metric_laplacian.eps'];
    saveas(gcf,fileName,'epsc');

    figure;
    set(gca,'FontSize',ts)
    plot(x,metric_histogram,'LineWidth',lw);
    xlabel('Image Shift (pixels)','FontSize',fs)
    ylabel('Histogram Energy','FontSize',fs)
    set(gca,'XTick',linspace(-5,5,3))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),ys))
    set(gca,'XMinorTick','on','YMinorTick','on')
    fileName = [char(pwd) '/figures/' 'metric_histogram.eps'];
    saveas(gcf,fileName,'epsc');

    figure;
    set(gca,'FontSize',ts)
    plot(x,metric_variance,'LineWidth',lw); 
    xlabel('Image Shift (pixels)','FontSize',fs)
    ylabel('Variance','FontSize',fs)
    set(gca,'XTick',linspace(-5,5,3))
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),ys))
    set(gca,'XMinorTick','on','YMinorTick','on')
    fileName = [char(pwd) '/figures/' 'metric_variance.eps'];
    saveas(gcf,fileName,'epsc');
  end
  
  %% compare to inverse radon transform - no translation, with translation, with rotation and translation



end