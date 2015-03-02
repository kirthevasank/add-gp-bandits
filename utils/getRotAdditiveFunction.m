function [func, funcProps] = getRotAdditiveFunction(numDims, numDimsPerGroup)

  % First get the additive Function
  [f, fProps] = getAdditiveFunction(numDims, numDimsPerGroup);

  % Now create an arbitrary rotation
  A = randn(numDims, numDims); A = orth(A);
  func = @(X) f(X * A);

  % Now copy the properties over
  funcProps.maxPt = fProps.maxPt * A;
  funcProps.maxVal = fProps.maxVal;
  funcProps.decomposition = fProps.decomposition;
  funcProps.A = A;

end

