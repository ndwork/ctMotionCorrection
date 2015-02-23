function transImg = translateImg( img, trans )
  % Function for translating an image in 2D
  % trans = [transY transX]

  [M, N] = size( img );
  
  hTrans = mod(trans(2),N);
  vTrans = mod(trans(1),M);

  % find x coordinates of translated image
  x = 1:N;
  transX = x - hTrans;
  % check for pixels falling off from right
  tmpX = transX > N;
  transX = transX - tmpX*N;
  % check for pixels falling off from left
  tmpX = transX < 1;
  transX = transX + tmpX*N;

  % find y coordinates of translated image
  y = 1:M;
  transY = y - vTrans;
  % check for pixels falling off from top
  tmpY = transY > M;
  transY = transY - tmpY*M;
  % check for pixels falling off from bottom
  tmpY = transY < 1;
  transY = transY + tmpY*M;
  
  padImgX = [img img(:,1)];

  x0 = floor(transX);
  x1 = ceil(transX);
  tmpX = repmat( transX-x0, [M 1] );
  transImgX = padImgX(:,x0) + (padImgX(:,x1)-padImgX(:,x0)) .* tmpX;
  
  padImgY = [ transImgX; transImgX(1,:); ];
  
  y0 = floor(transY);
  y1 = ceil(transY);
  tmpY = repmat( (transY-y0)', [1 N] );
  transImg = padImgY(y0,:) + (padImgY(y1,:)-padImgY(y0,:)) .* tmpY;

end
