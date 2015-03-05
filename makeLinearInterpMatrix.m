
function M = makeLinearInterpMatrix( x, y, xq )
  nq = numel(xq);
  ny = numel(y);
  M = zeros( nq, ny );

  minX = min(x);
  maxX = max(x);

  for i=1:nq
    if xq(i) < minX || xq(i) > maxX, continue; end;
    
    diff = xq(i) - x;
    [~,minDiffIndx] = min(diff);
    
    if xq(i) > x(minDiffIndx)
      b = xq(i) - x(minDiffIndx);
      a = x(minDiffIndx+1) - xq(i);
      xDist = x(minDiffIndx+1) - x(minDiffIndx);
      b = b / xDist;
      a = a / xDist;
      M(i,minDiffIndx) = a;
      M(i,minDiffIndx+1) = b;
    else
      b = xq(i) - x(minDiffIndx-1);
      a = x(minDiffIndx) - xq(i);
      xDist = x(minDiffIndx) - x(minDiffIndx-1);
      b = b / xDist;
      a = a / xDist;
      M(i,minDiffIndx-1) = a;
      M(i,minDiffIndx) = b;
    end

  end
end

