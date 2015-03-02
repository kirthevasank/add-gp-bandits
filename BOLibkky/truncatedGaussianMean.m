function truncMean = truncatedGaussianMean(mu, sigma, trunc)
% computes the value E[max(0, x)] where x~N(mu, sigma^2)

  y = mu - trunc;
  varZeroIdxs = (sigma == 0);

  truncMean = varZeroIdxs .* max(y, 0) + (~varZeroIdxs) .*  ...
    (y .* normcdf(y./sigma) + sigma .* normpdf(y./sigma) );

%   if sigma == 0
%     truncMean = max(y, 0);
%   else
%     truncMean = y * normcdf(y/sigma) + sigma * normpdf(y/sigma);
%   end
% 
end

