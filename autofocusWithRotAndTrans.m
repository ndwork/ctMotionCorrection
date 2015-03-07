
function [finaltranslation] = autofocusWithRotAndTrans(sinogram, ...
  minX, maxX, dx, minY, maxY, dy, maxRot, dr, nRows, nCols, nDetectors, ...
  pixSize, detSize, thetas, method)

  nThetas = numel(thetas);
  translations = zeros( nThetas, 2 );

  x = minX:dx:maxX; 
  y = minY:dy:maxY; 
  r = 0:dr:maxRot;

  metric= zeros(numel(r),numel(x),numel(y));

  for k = 1:numel(r)
    maxRotation = r(k) * pi/180;
    rotations = linspace(0,maxRotation,nThetas);

    for i = 1:numel(x)
      maxHorizontalShift = x(i);

      for j = 1:numel(y)
        maxVerticalShift = y(j);

        translations(:,1) = linspace(0,maxVerticalShift,nThetas);
        translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

        recon = ctCorrectForRotAndTrans( sinogram, nDetectors, ...
          detSize, thetas, rotations, translations, nCols, nRows, ...
          pixSize, 'method', method );

        metric(k,i,j) = normgrad(recon);
      end
    end
  end

  [rot,row,col] = find(metric==max(metric(:)));
  finaltranslation = [r(rot),y(col),x(row)];    
end

function [out] = normgrad(img)
  imgv = img(:);
  out = sum((abs(conv([1;-1],imgv))./sum(abs(conv([1;-1],imgv)))).^2);    
end
