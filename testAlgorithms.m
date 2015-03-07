function [recon,costs] = testAlgorithms( nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize, varargin)
  % Optional Variables:
  % method:
	%   'GD' for Gradient Descent
  %   'LADMM' for Linearized ADMM
  %   'PC' for Pock Chambolle

  defaultMethod = 'PC';
  defaultInput = 'rand';  
  p = inputParser;
  p.addOptional( 'method', defaultMethod );
  p.addOptional('inputMatrix', defaultInput);
  p.parse( varargin{:} );
  method = p.Results.method;
  inputMatrix = p.Results.inputMatrix;
  
  maxVerticalShift = max(translations(:,1));
  maxHorizontalShift = max(translations(:,2));
  switch inputMatrix
    
    case 'rand'
      random = rand(100,nRows);
      xTest = rand(nRows,nCols);
      xPad = padImgForRadon( xTest, maxHorizontalShift, maxVerticalShift, ...
        pixSize );

      applyE = @(u) random*u;
      applyET = @(u) random'*u;

      sinogram = applyE(xPad);
      
      gamma = 0;

    case 'deblur'
      % kernel = fspecial('gaussian',[15 15],7);
      kernel = fspecial('gaussian',[9 9],4); % use this one with the Barbara image below.
      % kernel = randn(15,15);

      kernelTrans = fliplr(flipud(kernel));

      applyE = @(x) imfilter(x,kernel,'circular');  % K is a blur operator.
      applyET = @(y) imfilter(y,kernelTrans,'circular');

      I = imread('barbara.png');
      I = double(I(:,:,1));
      I = imresize(I,.5);

      [nRows,nCols] = size(I);

      mn = min(I(:));
      I = I - mn;
      mx = max(I(:));
      I = I/mx;
      xTest = I;

      b = applyE(I);
      noise = (1e-3)*randn(nRows,nCols);
      b = b + noise;
      figure('Name','b')
      imshow(b,[]) % b is a blurry image, ready to be deblurred.
      sinogram = b;

      gamma = 1/300000; % This is the regularization parameter that gets multiplied by the regularization term TV(x).

  end
  
  
  switch method

    case 'GD'     % Gradient Descent
      costs = [];
      recon = ctCorrectForTranslation_GD_test( sinogram, ...
        applyE, applyET, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize );

    case 'LADMM'  % Linearized ADMM
      [recon, costs] = ctCorrectForTranslation_LADMM_test( sinogram, ...
        applyE, applyET, ...
        nDetectors, detSize, thetas, translations, nCols, nRows, ...
        pixSize, gamma );
      
      error = max(norm(recon(:) - xTest(:),2) / norm(xTest(:),2));
      disp(['Error is: ', num2str(error)]);
      disp(['gamma was set to: ', num2str(gamma)])

    case 'PC'     % Pock-Chambolle
      
      [recon, costs] = ctCorrectForTranslation_PC_test( sinogram, ...
        applyE, applyET, nDetectors, detSize, thetas, translations,...
        nCols, nRows, pixSize, gamma );
      
      figure; imshow( xTest, [] );
      title('Original Image')
      
      error = max( norm(recon(:) - xTest(:),2) / norm(xTest(:),2));
      disp(['Error is: ', num2str(error)]);
      disp(['gamma was set to: ', num2str(gamma)])
      
      

  end

end

