
function run_ctCorrectForTranslation
  clear; close all;
  addpath(genpath('.'));

  im = phantom();
  delta = 0.001;

  figure( 'name', 'Original Image' );
  imshow( im, [] );


  nDetectors = 500;
  dSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;

  % Reconstruction parameters
  cx = 0;   Nx=256;   dx=delta;
  cy = 0;   Ny=256;   dy=delta;


  type = 'fast';
  sinogram = ctRadon( im, delta, nDetectors, dSize, thetas, type );


end

