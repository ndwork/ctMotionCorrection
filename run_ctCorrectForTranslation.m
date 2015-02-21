
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
  nThetas = numel(thetas);

  maxVerticalShift = 0.01;
  maxHorizontalShift = 0.02;
  translations = zeros( nThetas, 2 );
%   translations(:,1) = linspace(0,maxVerticalShift,nThetas);
%   translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

  % Reconstruction parameters
  cx = 0;   Nx=256;   dx=delta;
  cy = 0;   Ny=256;   dy=delta;

  %sinogram = radonWithTranslation( im, delta, nDetectors, dSize, ...
  %  thetas, translations );
load 'junk.mat'

  recon = ctCorrectForTranslation( sinogram, delta, nDetectors, dSize, ...
    thetas, translations );


end

