
function out = RTWithTT( sinogram, transMatricesT, nCols, RT )
  %nCols = number of columns of the reconstructed image

  [N,~] = size(RT);
  nRows = N / nCols;
  nTranslations = numel(transMatricesT);


  out = zeros(nRows*nCols,1);
  parfor i=1:nTranslations
    sinoVec = zeros( size(sinogram) );
    sinoVec(i,:) = sinogram(i,:);
    sinoVec = sinoVec(:);

    bp = RT * sinoVec;   % back projected
    
    translated = transMatricesT{i} * bp;

    out = out + translated;
  end

  out = reshape( out, [nRows nCols] );
end

