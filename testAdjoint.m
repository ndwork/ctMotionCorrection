function [out, error] = testAdjoint(applyA,applyATrans,Mx,Nx,My,Ny)
  % function to test if applyATrans is the adjoint of applyA
  % Mx and Nx - size of test samples x
  % My and Ny - size of test samples y

  x = rand(Mx,Nx);
  y = rand(My,Ny);

  Ax = applyA(x);
  ATy = applyATrans(y);

  prod1 = sum(sum(Ax.*y));
  prod2 = sum(sum(x.*ATy));
  
  error = abs(prod1 - prod2);
  if error < 1e-10
    out = 1;
  else
    out = 0;
  end
end