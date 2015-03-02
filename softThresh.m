
function out = softThresh( in, thresh )

  tmp = abs(in) - thresh;
  out = sign(in) .* max( tmp, 0 );

end