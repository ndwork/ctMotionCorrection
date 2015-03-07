
function run_testAlgorithms
  clear; close all;
  addpath(genpath('.'));

  % Reconstruction parameters
  method = 'PC';    % Options: GD, PC, LADMM
  inputMatrix = 'deblur'; % Options: rand (random matrix), deblur (image 
                        % deblurring function
  cx = 0;   nCols=50;
  cy = 0;   nRows=30;
  pixSize = 0.001; % meters / pixel

  detSize = 0.001;
  dTheta = 2 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  nThetas = numel(thetas);

  nDetectors = nCols*2;

%   maxVerticalShift = 0.01; % in meters
%   maxHorizontalShift = 0.02; % in meters

maxVerticalShift = 0; % in meters
  maxHorizontalShift = 0; % in meters
  translations = zeros( nThetas, 2 );
  %translations(:,1) = linspace(0,maxVerticalShift,nThetas);
  %translations(:,2) = linspace(0,maxHorizontalShift,nThetas);

  %profile on
  tic;
  [recon,costs] = testAlgorithms( nDetectors, detSize, ...
    thetas, translations, nCols, nRows, pixSize, 'method', method, ...
    'inputMatrix',inputMatrix);
  timeTaken = toc;
  %profile off
  %profile viewer

  disp(['Time taken: ', num2str(timeTaken)]);
%   figure; imshow( imresize(recon,10,'nearest'), [] );
  figure; imshow( recon, [] );
  title('Reconstructed image');
  
  figure; plot( costs, 'LineWidth', 2 );
  xlabel('Iteration'); ylabel('Cost Function');

end

