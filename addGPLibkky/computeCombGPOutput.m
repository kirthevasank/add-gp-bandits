function [yMu, yStd, yK] = computeCombGPOutput(Xte, Xtr, decomposition, L, alpha, ...
  bws, scales, meanFuncs, commonMeanFunc)
% Returns the predictive mean and variance for the combined GP

  numGroups = numel(decomposition);
  K12 = combinedKernel(Xtr, Xte, decomposition, bws, scales);
  K22 = combinedKernel(Xte, Xte, decomposition, bws, scales);

  % Compute the outputs
  yMu = combinedMeanFunc(Xte, commonMeanFunc, meanFuncs, decomposition) + ...
        K12' * alpha;
  V = L \ K12;
  yK = K22 - V'*V;
  yStd = sqrt(real(diag(yK)));
end
