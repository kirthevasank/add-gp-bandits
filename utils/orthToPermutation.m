function [P, Z] = orthToPermutation(A, costIdx, pwr)
% Takes in an Orthogonal matrix and returns the "closest" permutation matrix.

  if ~exist('costIdx', 'var') | isempty(costIdx)
    costIdx = 2; % Since it corresponds to a norm
    % for instance if pwr =1, l1 norm if pwr=2, frobenius norm etc.
  end

  if ~exist('pwr', 'var') | isempty(pwr)
    pwr = 2;
  end

  % prelims
  D = size(A, 1);
  p = size(A, 2);
  if D ~= p
    [A, ~] = qr(A);
  end

  % First get the Cost matrix for the assignment
  C = getHungarianCostMatrixOrth(A, costIdx, pwr);

  % Now use Hungarian
  Z = hungarian(C);
  Z = Z(1:p);

  % Now construct the permutation matrix
  P = zeros(D, p);
  for i = 1:p
    P(Z(i), i) = 1;
  end

end

