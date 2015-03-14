
function sinogram = RWithTranslation( img, translations_pix, nDet, RT )
  % sinogram = RWithTranslation( img, translations_pix, nDet, R )
  % nDet is the number of detectors

  %[M,~] = size(R);
  [~,M] = size(RT);
  nThetas = M / nDet;
  sinogram = zeros(nThetas,nDet);

  [nTranslations,~] = size(translations_pix);
  if nTranslations ~= nThetas, error('Input error'); end;

  parfor th=1:nThetas
    thisTrans_pix = translations_pix(th,:);
    translated = translateImg( img, thisTrans_pix );

    %tmpSino = R * translated(:);
    %tmpSino = reshape( tmpSino, [nThetas nDet] );
    %sinogram(th,:) = tmpSino(th,:);

    indxs1D = th + ((1:nDet)-1)*nThetas;
    results = RT(:,indxs1D)' * translated(:);
    %results = R(indxs1D,:) * translated(:);
    sinogram(th,:) = results;
  end

end
