function [out] = testAdjoint(applyE,applyETrans,Mx,My,Nx,Ny)
% function to test if applyETrans applies the adjoint of applyE
% Mx and Nx - size of test samples x
% My and Ny - size of test samples y

  x = rand(Mx,Nx);
  y = rand(My,Ny);

  tmp1 = dot(applyE(x),y);
  tmp2 = dot(x,applyETrans(y));
  
  if abs(tmp1 - tmp2) < 1e-10
    out = 1;
  else
    out = 0;
  end


end