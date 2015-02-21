
function run_ctCorrectForTranslation
  clear; close all;
  addpath(genpath('.'));

  im = phantom();

  nDetectors = 500;
  detSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);
  
  % Reconstruction parameters
  cx = 0;   Nx=256;
  cy = 0;   Ny=256;
  pixSize = 0.001; % meters / pixel

  maxVerticalShift = 0.01; % in meters
  maxHorizontalShift = 0.02; % in meters
  translations = zeros( nThetas, 2 );
  translations(:,1) = linspace(0,maxVerticalShift,nThetas);
  translations(:,2) = linspace(0,maxHorizontalShift,nThetas);
  
  % pad the phantom image so that it always stays in the field of view
  xShiftPix = maxHorizontalShift / pixSize;
  yShiftPix = maxVerticalShift / pixSize;
  if xShiftPix >= 0 
    im = [im zeros(size(im,1),ceil(xShiftPix))];
  else
    im = [zeros(size(im,1),ceil(abs(xShiftPix))) im];
  end
  if yShiftPix >= 0
    im = [zeros(ceil(yShiftPix),size(im,2)); im];
  else
    im = [im; zeros(ceil(abs(yShiftPix)),size(im,2))];
  end

  %sinogram = radonWithTranslation( im, pixSize, nDetectors, detSize, ...
  %  thetas, translations );
load 'phSinogram.mat'

  recon = ctCorrectForTranslation( sinogram, nDetectors, detSize, ...
    thetas, translations, Nx, Ny, pixSize );


end



