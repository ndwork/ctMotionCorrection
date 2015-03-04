
function recon = ctCorrectForTranslation_GD( sinogram, nDetectors, ...
  detSize, thetas, translations, nCols, nRows, pixSize )
  % This function uses Pock-Chambolle to determine the reconstruction
  % image based on the known translations
  % sinogram is an MxN array
  % translation is an Mx2 element array; each row of the array is the
  %   translation for the corresponding row of the sinogram

  if nargin < 10, lambda=1; end;

  %gamma = 1d-5;   % Regularization parameter
  gamma = 0;

  applyD1 = @(u) cat(2, u(:,2:end) - u(:,1:end-1), zeros(nRows,1));
  applyD2 = @(u) cat(1, u(2:end,:) - u(1:end-1,:), zeros(1,nCols));
  applyD1T = @(u) cat(2, -u(:,1), u(:,1:end-2) - u(:,2:end-1), u(:,end-1));
  applyD2T = @(u) cat(1, -u(1,:), u(1:end-2,:) - u(2:end-1,:), u(end-1,:));

  applyE = @(u) radonWithTranslation( u, pixSize, nDetectors, ...
    detSize, thetas, translations );

  cx = 0;  cy = 0;
  applyET = @(u) backprojectionWithTranslation( u, thetas, ...
    detSize, cx, cy, nCols, nRows, pixSize, translations );

  recon = lsqr( @applyForGD, sinogram(:) );

end

function out= applyForGD( in, type )
  switch type
    case 'transp'
      out = applyKT( in{1}, in{2}, in{3} );
    case 'notransp'
      out = applyK( in );
  end
end


function out = applyK( uE, uD1, uD2 )
  out = cell(3,1);
  out{1} = applyE( uE );
  out{2} = applyD1( uD1 );
  out{3} = applyD2( uD2 );
end

function out = applyKT( uE, uD1, uD2 )
  outE = applyET( uE );
  outD1 = applyD1T( uD1 );
  outD2 = applyD2T( uD2 );
  out = outE + outD1 + outD2;
end

