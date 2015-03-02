function plotAddGPOutputs(f, funcHs, trueF, trueFuncHs, decomposition, bounds)
% This is just a general utility function to plot the outputs of an additive GP.

  numGroups = numel(decomposition);
  numDims = size(bounds, 1);

  for k = 1:numGroups

    coords = decomposition{k};
    currBounds = bounds(coords, :);
    if numel(coords) == 1
      figure;
      t = linspace(currBounds(1), currBounds(2), 200)';
      [mu, sigma] = funcHs{k}(t);
      plot(t, mu, 'b'); hold on,
      plot(t, trueFuncHs{k}(t), 'r'); hold on,
      plot(t, mu + sigma, 'g');
      plot(t, mu - sigma, 'g');
      legend('Estimated', 'True');
      titlestr = sprintf('Func-Idx: %d', k);
      title(titlestr);
    elseif numel(coords) == 2
      figure;
      plot2DFunction(funcHs{k}, currBounds, 'surfc', 'b'); hold on,
      plot2DFunction(trueFuncHs{k}, currBounds, 'surfc', 'r');
      legend('Estimated', 'True');
      titlestr = sprintf('Func-Idx: %d', k);
      title(titlestr);
    end
  end

  % Now plot the complete function
  if numDims == 2
    figure;
    meanF = @(t) getMean(t, f);
    plot2DFunction(meanF, bounds, 'contour', 'b'); hold on,
    plot2DFunction(trueF, bounds, 'contour', 'r');
    legend('Estimated', 'True');
    title('Complete Function');
  end
end

function mu = getMean(t, f)
  [mu, ~] = f(t);
end
