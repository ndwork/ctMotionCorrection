
function [recon] = filteredBackProjection(sinogram, thetas, detSize, ...
  nCols, nRows, pixSize)

  nThetas = numel(thetas);
  nDetectors = size(sinogram,2);

  dLocIdx = [-(nDetectors-1)/2 : 1 : (nDetectors-1) / 2];
  
  fftSino = fftshift(fft(ifftshift(sinogram,2),[],2),2);

  dk = 1 / numel(dLocIdx);
  N = nThetas;
  
  rhoFilter = ((((dk^2)*pi) / N) * abs(dLocIdx));
  rhoFilter( dLocIdx==0 ) = 0.25 * ((dk^2 * pi) / N);
  filter = ones(nThetas,1) * rhoFilter;
  
  newFFTSino = filter .* fftSino;
  
  newSino = fftshift(ifft(ifftshift(newFFTSino,2),[],2),2);
  newSino = real(newSino);
  
  recon = backprojection(newSino, -thetas, detSize, ...
    nCols, nRows, pixSize);

end
