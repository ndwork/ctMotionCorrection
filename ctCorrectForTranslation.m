
function recon = ctCorrectFortranslation( sinogram, delta, nDetectors, ...
  dSize, thetas, translations )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  x = zeros( size(sinogram) );
  y = zeros( size(sinogram) );
  xbar = zeros( size(sinogram) );

  nIter = 1000;
  
  for i=1:nIter
    y = 
  end
  
  nrmK = estimateNormByPowerIteration(appK,appKTrans,randn(M,1),1000);
  sigma = 1/(tau*nrmK^2);
  
  
end




