function [func, funcProperties] = get2Modal1DFunction
% Returns a function of 2 modes in 1D

  FUNC_LB = -700;
  bounds = [-1, 1];
  funcProperties.bounds = bounds;
  gaussVar = 0.01;
  funcProperties.gaussVar = gaussVar;
  centres = [-0.5; 0.4];
  funcProperties.centres = centres;

  probs = [0.15; 0.85];
  funcProperties.centreProbs = probs;
  func = @(t) max( - 700, log( ...
    probs(1) * mvnpdf(t, centres(1), gaussVar) + ...
    probs(2) * mvnpdf(t, centres(2), gaussVar) ) );
  [~,maxCentreIdx] = max(probs);
  funcProperties.maxPt = funcProperties.centres(maxCentreIdx, :);
  funcProperties.maxVal = func(funcProperties.maxPt);

  if funcProperties.maxVal < FUNC_LB
    error('Decrease Lower Bound. MaxVal = %.4f\n', funcProperties.maxVal);
  end

end

