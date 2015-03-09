
function sinogram = RWithT( img, transMatrices, nDet, R )
  % nDet is the number of detectors

  [M,~] = size(R);
  nThetas = M / nDet;
  sinogram = zeros(nThetas,nDet);

  nTranslations = numel( transMatrices );
  if nTranslations ~= nThetas, error('Input error'); end;
  
  for i=1:nThetas
    translated = transMatrices{i} * img(:);
    tmpSino = R * translated;
    tmpSino = reshape( tmpSino, [nThetas nDet] );
    sinogram(i,:) = tmpSino(i,:);
  end

end

