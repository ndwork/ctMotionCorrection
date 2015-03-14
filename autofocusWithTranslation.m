<<<<<<< HEAD

function [bestTranslations,metric] = autofocusWithTranslation(sinogram, maxX, dx, ...
  maxY, dy, nRows, nCols, pixSize, detSize, thetas, method)

  [nThetas,nDetectors] = size(sinogram);

  x = [-maxX:dx:maxX]; 
  y = [-maxY:dy:maxY]; 

  metric = zeros(numel(x),numel(y));
  translations = zeros( nThetas, 2 );

  %R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
  %  detSize, thetas);
%save( 'RadonMatrix_autofocus.mat', 'R' );
load 'RadonMatrix_autofocus.mat';

  nX = numel(x);
  for i = 1:nX
    disp(['Working on ', num2str(i), ' of ', num2str(nX)]);

    maxHorizontalShift = x(i); % in meters
    translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

    for j = 1:numel(y)
      maxVerticalShift = y(j); % in meters
      translations(:,1) = linspace(0,maxVerticalShift,nThetas);

      [recon,~] = ctCorrectForTranslation( sinogram, nDetectors, ...
        detSize, thetas, translations, nCols, nRows, pixSize, ...
        'method', method, 'radonMatrix', R );

      metric(i,j) = normgrad(recon);
    end
  end

  [row,col] = find( metric==max(metric(:)) );
  bestTranslations = [y(col),x(row)];  
end

function [out] = normgrad(img)
  imgv = img(:);
  out = sum((abs(conv2([1;-1],imgv))./sum(abs(conv2([1;-1],imgv)))).^2);    
end
=======
function [final] = autofocusWithTranslation(sinogram, minX, maxX, dx, minY, maxY, dy, nRows, nCols, nDetectors, pixSize, detSize, thetas, method)
    nThetas = numel(thetas);

    x = [minX:dx:maxX]; 
    y = [minY:dy:maxY]; 
        
    metric = zeros(numel(x),numel(y));
    for i = 1:numel(x) 
        maxHorizontalShift = x(i); % in meters
        for j = 1:numel(y)          
            maxVerticalShift = y(j); % in meters
            translations = zeros( nThetas, 2 );
            translations(:,1) = linspace(0,maxVerticalShift,nThetas);
            translations(:,2) = linspace(0,maxHorizontalShift,nThetas);
        
            [recon,costs] = ctCorrectForTranslation( sinogram, nDetectors, detSize, ...
            thetas, translations, nCols, nRows, pixSize, 'method', method );

            metric(i,j) = normgrad(recon);
        end
    end
    [row,col] = find(metric==max(metric(:)));
    final = [y(col),x(row)];  
end
function [out] = normgrad(img)
    imgv = img(:);
    out = sum((abs(conv([1;-1],imgv))./sum(abs(conv([1;-1],imgv)))).^2);    
end
>>>>>>> 0e77338396d54c693abca30d1889dfcb39347bbb
