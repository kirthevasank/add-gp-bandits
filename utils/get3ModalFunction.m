function [func, funcProperties] = get3ModalFunction(numDims)
% Returns a function handle with 3 modes in numDims dimensions. The function has
% 3 modes (see below) and the maximum occurs at the following
% numDims-dimensional vec: [0.59; 0.79; (numDims-2) times 0.41]

  FUNC_LB = -700;

  bounds = repmat([-1, 1], numDims, 1);
  funcProperties.bounds = bounds;
  funcProperties.gaussVar = 0.01 * numDims^0.1;
  funcProperties.covarDiag = funcProperties.gaussVar * ones(1, numDims);
  funcProperties.centres12 = [0.62 -0.38; 0.18 0.58; -0.58 -0.56];
  funcProperties.centresRest = [-0.66 -0.19 0.62]';
  funcProperties.centres = ...
    [funcProperties.centres12, repmat(funcProperties.centresRest, 1, numDims -2)];
  probs = [0.1; 0.8; 0.1];
%   probs = [0.3; 0.5; 0.3];
  funcProperties.centreProbs = probs;
  func = @(t) max( FUNC_LB, log( ...
    probs(1) * mvnpdf(t, funcProperties.centres(1, :), funcProperties.covarDiag ) +...
    probs(2) * mvnpdf(t, funcProperties.centres(2, :), funcProperties.covarDiag ) +...
    probs(3) * mvnpdf(t, funcProperties.centres(3, :), funcProperties.covarDiag ) ) ...
    );
  [~, maxCentreIdx] = max(probs);
  funcProperties.maxPt = funcProperties.centres(maxCentreIdx, :);
  funcProperties.maxVal = func(funcProperties.maxPt);

  if funcProperties.maxVal < FUNC_LB
    error('Decrease Lower Bound. MaxVal = %.4f\n', funcProperties.maxVal);
  end

end

