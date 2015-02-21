
function run_ctCorrectForTranslation
  clear; close all;
  addpath(genpath('.'));

  im = phantom();

  nDetectors = 500;
  detSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);

  maxVerticalShift = 0.01;
  maxHorizontalShift = 0.02;
  translations = zeros( nThetas, 2 );
%   translations(:,1) = linspace(0,maxVerticalShift,nThetas);
%   translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

  % Reconstruction parameters
  cx = 0;   Nx=256;
  cy = 0;   Ny=256;
  pixSize = 0.001;

  %sinogram = radonWithTranslation( im, pixSize, nDetectors, detSize, ...
  %  thetas, translations );
load 'phSinogram.mat'

  recon = ctCorrectForTranslation( sinogram, nDetectors, detSize, ...
    thetas, translations, Nx, Ny, pixSize );


end

