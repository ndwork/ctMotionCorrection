
function [bestTranslations,metric] = autofocusWithRotAndTrans(sinogram, ...
  maxX, dx, maxY, dy, maxRot, dr, nRows, nCols, pixSize, detSize, ...
  thetas, method)

  [nThetas,nDetectors] = size(sinogram);

  x = [-maxX:dx:maxX]; 
  y = [-maxY:dy:maxY]; 

  metric = zeros(numel(x),numel(y));
  translations = zeros( nThetas, 2 );

  R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
    detSize, thetas);

  nR = numel(r);
  for k = 1:nR
    maxRotation = r(k) * pi/180;
    rotations = linspace(0,maxRotation,nThetas);

    nX = numel(x);
    for i = 1:nX
      disp(['Working on ', num2str(i), ' of ', num2str(nX)]);

      maxHorizontalShift = x(i); % in meters
      translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

      for j = 1:numel(y)
        maxVerticalShift = y(j); % in meters
        translations(:,1) = linspace(0,maxVerticalShift,nThetas);

        recon = ctCorrectForRotAndTrans( sinogram, nDetectors, ...
          detSize, thetas, rotations, translations, nCols, nRows, ...
          pixSize, 'method', method, 'radonMatrix', R );

        metric(i,j) = normgrad(recon);
      end
    end
  end

  [row,col] = find( metric==max(metric(:)) );
  bestTranslations = [y(col),x(row)];  
end

function [out] = normgrad(img)
  imgv = img(:);
  out = sum((abs(conv2([1;-1],imgv))./sum(abs(conv2([1;-1],imgv)))).^2);    
end
