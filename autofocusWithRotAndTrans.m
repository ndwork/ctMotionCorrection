function [finaltranslation] = autofocusWithRotAndTrans(sinogram, minX, maxX, dx, minY, maxY, dy, maxRot, dr, nRows, nCols, nDetectors, pixSize, detSize, thetas, method)
    nThetas = numel(thetas);

    x = [minX:dx:maxX]; 
    y = [minY:dy:maxY]; 
    r = [0:dr:maxRot];
        
    metric= zeros(numel(x),numel(y),numel(r));
    
    for i = 1:numel(x)
        for j = 1:numel(y)
            for k = 1:numel(r)
                maxHorizontalShift = x(i);
                maxVerticalShift = y(j);
                translations = zeros( nThetas, 2 );
                translations(:,1) = linspace(0,maxVerticalShift,nThetas);
                translations(:,2) = linspace(0,maxHorizontalShift,nThetas);
                
                maxRotation = r(k) * pi/180;
                rotations = linspace(0,maxRotation,nThetas);
                
                recon = ctCorrectForRotAndTrans( sinogram, nDetectors, ...
       detSize, thetas, rotations, translations, nCols, nRows, pixSize, 'method', method );

                metric(i,j,k) = normgrad(recon);
            end
        end
    end
    [row,col,rot] = find(metric==max(metric(:)));
    finaltranslation = [r(rot),y(col),x(row)];    
end

function [out] = normgrad(img)
    imgv = img(:);
    out = sum((abs(conv([1;-1],imgv))./sum(abs(conv([1;-1],imgv)))).^2);    
end