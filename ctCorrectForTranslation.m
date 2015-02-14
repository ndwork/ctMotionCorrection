
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
%     y = 
  end
  
  K = [3 0 0;
       0 2 0;
       0 0 1];
  
  M = size(K,2);
  
  applyK = @(x) K*x;
  applyKTrans = @(x) K'*x;
  
  nrmK = estimateNormByPowerIteration(applyK,applyKTrans,randn(M,1),1000);
  sigma = 1/(tau*nrmK^2);
  
  
end




