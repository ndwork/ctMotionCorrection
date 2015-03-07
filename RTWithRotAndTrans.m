
function bp = RTWithRotAndTrans( sinogram, rotations, translations, ...
  nCols, RT )
  %nCols = number of columns of the reconstructed image

  [N,~] = size(RT);
  nRows = N / nCols;
  [nTranslations,~] = size(translations);


  bp = zeros(nRows,nCols);
  parfor i=1:nTranslations
    sinoVec = zeros( size(sinogram) );
    sinoVec(i,:) = sinogram(i,:);
    sinoVec = sinoVec(:);

    tmp = RT * sinoVec;   % back projected
    tmp = reshape( tmp, [nRows nCols] );

    thisTrans = translations(i,:);
    translated = translateImg( tmp, -thisTrans );

    thisRot = rotations(i);
    rotated = imrotate( translated, -thisRot, 'bilinear', 'crop' );
    
    bp = bp + rotated;
  end

end
