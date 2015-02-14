
function recon = ctCorrectFortranslation( sinogram, delta, nDetectors, ...
  dSize, thetas, translations )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1,:), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1Trans = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2Trans = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));


  sizeSino = size(sinogram);
  xE = zeros( sizeSino(2) );
  xD1 = zeros( sizeSino(2) );
  xD2 = zeros( sizeSino(2) );
  xBarE = zeros( sizeSino );
  xBarD1 = zeros( sizeSino );
  xBarD2 = zeros( sizeSino );
  yE = zeros( sizeSino );
  yD1 = zeros( sizeSino(2) );
  yD2 = zeros( sizeSino(2) );
  

  nIter = 1000;
  
  for i=1:nIter

  end
  
  
  
end




