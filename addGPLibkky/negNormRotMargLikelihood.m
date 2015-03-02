function [negNlml, nAG] = negNormMargLikelihood(sigmaSms, sigmaPrs, ...
  decomposition, A, X, y, meanFuncs, commonMeanFunc, noises, commonNoise) 

  nlml = normRotMargLikelihood(sigmaSms, sigmaPrs, decomposition, A, ...
    X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
  negNlml = -nlml;

end

