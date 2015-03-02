function util = getUCBUtility(x, funcH, numEvals)

  % Prelims
  numDims = size(x, 2); % Expecting each x to be a row vector here.

  % Set beta_t. Using Recommendation from Section 6 in Srinivas et al, ICML 2010
  delta = 0.01; % something we need to set.
  t = numEvals + 1;

%   % Simple UCB Rule
%   beta_t = 2 * log( numDims * (t*pi)^2 / (6 * delta) ) / 5;

  % Linear in D, log in t
  beta_t = numDims * log( 2*t)/5;

%   % UCB Rule for finite D. Here we use |D| = 1000;
%   beta_t = 2 *(log(1000) + log( (t*pi)^2 / (6*delta) )) /5;

%   % UCB Rule for f with bounded RKHS norm
%   beta_t = 2 * 10 + 300 * log(t)^(numDims + 1) * (log(t/delta))^3;

  % UCB rule for f\sim GP
%   beta_t = 2 *log(t^2 * 2*pi^2/(3*delta)) + ...
%           2*numDims*log(t^2 *numDims * 2 * sqrt(4*numDims/delta) );

  % Obtain mean and standard deviation
  [mu, sigma] = funcH(x);
  util = mu + sqrt(beta_t) * sigma;

end

