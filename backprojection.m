function [backProj] = backprojection(sinogram, thetas, detSize, ...
  nCols, nRows, pixSize)
  
  thetas_deg = (180 / pi) * thetas;
  nThetas = numel(thetas);
  nDetectors = size(sinogram,2);
  nPix = nCols;

  dLocIdx = [-(nDetectors-1)/2 : 1 : (nDetectors-1) / 2];
  dLocs = dLocIdx * detSize;  
  pLocs = [-(nPix-1)/2 : 1 : (nPix-1) / 2] * pixSize;

  backProj = zeros( nRows, nCols);
  
  
  xs = ones(nRows,1) * (1:nCols);
  ys = (1:nRows)' * ones(1,nCols);
  halfY = nRows/2;  ys = ys - halfY;
  halfX = nCols/2;  xs = xs - halfX;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(nCols/2,nRows/2);

  for k = 1:nThetas;
    interped = interp1(dLocs, sinogram(k,:), pLocs) * pixSize;
    
    % smear row
    smear = ones(nRows,1)*interped;
    
    masked = smear .* radiusMask;
    
    rotated = imrotate(masked,thetas_deg(k),'bilinear','crop');
    
    
    
    backProj = backProj + rotated;
  end

end