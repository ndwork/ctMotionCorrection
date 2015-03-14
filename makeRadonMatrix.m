
function R = makeRadonMatrix( nCols, nRows, pixSize, nDetectors, ...
  detSize, thetas)

  nThetas = numel(thetas);
  R = sparse(nDetectors*nThetas,nCols*nRows);
  
  parfor i=1:nCols*nRows
    if mod(i,100)==0
      disp(['Making Radon Matrix ', num2str(i), ' of ', ...
        num2str(nRows*nCols)]);
    end

    img = zeros(nRows,nCols);
    img(i) = 1;
    Rimg = radonTransform( img, pixSize, nDetectors, ...
      detSize, thetas );
    R(:,i) = Rimg(:);
  end

end
