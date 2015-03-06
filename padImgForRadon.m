function out = padImgForRadon( im, maxHorizontalShift, ...
  maxVerticalShift, pixSize )
  % pad the image so that it always stays in the field of view

  xShiftPix = maxHorizontalShift / pixSize;
  xPadding = zeros(size(im,1),ceil(abs(xShiftPix)));
  out = [xPadding im xPadding];
  yShiftPix = maxVerticalShift / pixSize;
  yPadding = zeros(ceil(abs(yShiftPix)),size(out,2));
  out = [ yPadding; out; yPadding ];
end
