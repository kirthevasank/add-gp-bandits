function [yMu, yStd, yK] = computeIndGPOutput(Xte, Xtr, coords, L, alpha, ...
  bw, scale, meanFunc)
% This function returns the predictive mean and GP for an individual GP

  % Compute K12 and K22
  K12 = augKernel(Xtr, Xte, coords, bw, scale);
  K22 = augKernel(Xte, Xte, coords, bw, scale);
  % Note that we are NOT adding noise to K22 here. 
  % If we have already observed at Xte then we need to account for each
  % individual noise. But if we are interested in prediction at a *new
  % unobserved* point then we need not add noise.
  % TODO: Write a function for this.

  % Now compute the outputs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  yMu = meanFunc(Xte) + K12' * alpha; % Predictive Mean
  V = L \ K12;
  yK = K22 - V'*V; % Predictive Variance
  yStd = sqrt(real(diag(yK)));
end
