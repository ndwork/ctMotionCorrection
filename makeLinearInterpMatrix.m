
function M = makeLinearInterpMatrix( x, xq )
  nq = numel(xq);
  nx = numel(x);
  M = zeros( nq, nx );

  minX = min(x);
  maxX = max(x);

  for i=1:nq
    if xq(i) < minX || xq(i) > maxX, continue; end;

    diff = xq(i) - x;
    [~,minDiffIndx] = min(abs(diff));
    
    if xq(i) == x(minDiffIndx)
      M(i,minDiffIndx) = 1;
    elseif xq(i) > x(minDiffIndx)
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

