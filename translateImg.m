function [transImg] = translateImg(img,trans)
% Function for translating an image in 2D
% trans = [transY transX]
  
  numX = size(img,2);
  numY = size(img,1);

  hTrans = mod(trans(2),numX);
  vTrans = mod(trans(1),numY);
  
  % x and y coordinates of original image
  x = [1:1:numX];
  y = [1:1:numY];
  [x, y] = meshgrid(x,y);
  
  % pad the x coords, y coords, and data
  padX = [x (numX+1)*ones(size(x,1),1)];
  padX = [padX; padX(1,:)];
  
  padY = [y; (numY+1)*ones(1,size(y,2))];
  padY = [padY padY(:,1)];
  
  padImg = [img img(:,1)];
  padImg = [padImg; padImg(1,:)];
  
  % find x coordinates of translated image
  transX = x - hTrans;
  % check for pixels falling off from right
  tmpX = transX > numX;
  transX = transX - tmpX*numX;
  % check for pixels falling off from left
  tmpX = transX < 1;
  transX = transX + tmpX*numX;
  
  % find y coordinates of translated image
  transY = y - vTrans;
  % check for pixels falling off from top
  tmpY = transY > numY;
  transY = transY - tmpY*numY;
  % check for pixels falling off from bottom
  tmpY = transY < 1;
  transY = transY + tmpY*numY;
  
  % translate image
  transImg = griddata(padX,padY,padImg,transX,transY);


end