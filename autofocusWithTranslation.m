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
