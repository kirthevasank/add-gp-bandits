function [nlml] = normRotMargLikelihood(sigmaSms, sigmaPrs, decomposition, A, ...
  X, y, meanFuncs, commonMeanFunc, noises, commonNoise)
% Returns the normalized marginal likelihood.
% Decomposition is the decomposition after applying Z = X*A;

  % prelims
  numPts = size(X, 1);
  numGroups = numel(decomposition);
  D = size(A, 1);
  p = size(A, 2);

  % apply transformation and compute normalized marginal likelihood
  Z = X * A;
  Ky = combinedKernelNoise(Z, decomposition, sigmaSms, sigmaPrs, noises, ...
        commonNoise);
  L = stableCholesky(Ky);
  y_ = y - combinedRotMeanFunc(X, Z, commonMeanFunc, meanFuncs, decomposition);
  alpha = L' \ (L \ y);
  nlml = -1/2 * y_' * alpha - sum(log(diag(L))) - numPts/2 * log(2*pi);

end


function mu0 = combinedRotMeanFunc(X, Z, commonMeanFunc, meanFuncs, ...
                 decomposition)
% The common Mean Func takes in X as its arguments while meanFuncs take in Z =
% X*A as its arguments

  mu0 = commonMeanFunc(X);
  numGroups = numel(decomposition);
  for k = 1:numGroups
    coords = decomposition{k};
    mu0 = mu0 + meanFuncs{k}( Z(:, coords) );
  end

end
