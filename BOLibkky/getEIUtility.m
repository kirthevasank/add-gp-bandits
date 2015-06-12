function util = getEIUtility(x, gpFuncH, trunc)
% Expected Improvement Utility. Applies only to non-additive functions.
  [mu, s] = gpFuncH(x);
  util = truncatedGaussianMean(mu, s, trunc);
end


function truncMean = truncatedGaussianMean(mu, sigma, trunc)
% computes the value E[max(0, x)] where x~N(mu, sigma^2)

  y = mu - trunc;
  varZeroIdxs = (sigma == 0);

  truncMean = varZeroIdxs .* max(y, 0) + (~varZeroIdxs) .*  ...
    (y .* normcdf(y./sigma) + sigma .* normpdf(y./sigma) );

end

