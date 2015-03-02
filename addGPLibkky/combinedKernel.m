function K = combinedKernel(X1, X2, decomposition, bws, scales)
% This is the complete Kernel K0 = sum_i Ki
  numGroups = numel(decomposition);
  n1 = size(X1, 1);
  n2 = size(X2, 1);
  K = zeros(n1, n2);
  for k = 1:numGroups
    coords = decomposition{k};
    bw = bws(k);
    scale = scales(k);
    K = K + augKernel(X1, X2, coords, bw, scale);
  end
end
