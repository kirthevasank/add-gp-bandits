function mu0 = combinedMeanFunc(X, commonMeanFunc, meanFuncs, decomposition)
% The total mean function (obtained by adding the common and individual mean
% Functions)

  mu0 = commonMeanFunc(X);
  numGroups = numel(decomposition);
  for k = 1:numGroups
    coords = decomposition{k};
    mu0 = mu0 + meanFuncs{k}(X(:, coords) );
  end
end

