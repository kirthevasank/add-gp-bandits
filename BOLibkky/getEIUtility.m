function util = getEIUtility(x, gpFuncH, trunc)
% Expected Improvement Utility
  [mu, s] = gpFuncH(x);
  util = truncatedGaussianMean(mu, s, trunc);
end
