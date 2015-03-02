function K = GaussKernel(h, X, Y)
% Returns the Kernel Matrix for a Gaussian Kernel of bandwidth h.
% X is an nxd matrix. K is an nxn matrix.
% If Y is nonempty then returns the gaussian kernel for XxY

  if ~exist('Y', 'var') | isempty(Y)
    Y = X;
  end

  d = size(X, 2); % dimensions
  D = dist2(X, Y);
  K = 1/sqrt(2*pi * h^2)^(d) * exp(-D/(2*h^2));

end
