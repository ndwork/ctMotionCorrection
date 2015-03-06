
function sinogram = RWithTranslation( img, translations_pix, nDet, R )
  % nDet is the number of detectors

  [M,~] = size(R);
  nThetas = M / nDet;
  sinogram = zeros(nThetas,nDet);

  [nTranslations,~] = size(translations_pix);
  if nTranslations ~= nThetas, error('Input error'); end;

  for i=1:nThetas
    thisTrans_pix = translations_pix(i,:);
    translated = translateImg( img, thisTrans_pix );
    %sinogram(i,:) = R(nDet*(i-1)+1:nDet*i,:) * translated(:);
    tmpSino = R * translated(:);
    tmpSino = reshape( tmpSino, [nThetas nDet] );
    sinogram(i,:) = tmpSino(i,:);
  end

end
