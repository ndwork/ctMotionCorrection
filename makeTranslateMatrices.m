
function out = makeTranslateMatrices( M, N, trans )

  [nTrans,~] = size(trans);
  out = cell( nTrans, 1 );

  parfor k=1:nTrans
    if mod(k,10)==0
      disp(['Working on translation matrix ', num2str(k), ' of ', ...
        num2str(nTrans)]);
    end

    thisMatrix = sparse(M*N,M*N);

    for i=1:M*N
      basisVec = zeros(M,N);
      basisVec(i) = 1;
      translated = translateImg( basisVec, trans(k,:) );
      thisMatrix(:,k) = translated(:);
    end

    out{k} = thisMatrix;
  end

end
