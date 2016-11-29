function [mu, KPost, Mus, KPosts, combinedFuncH, funcHs, ...
  sigmaSmOpts, sigmaPrOpts] = ...
  addGPMargLikelihood(X, y, Xtest, decomposition, hyperParams)
% Pick the parameters for the GP by maximizing the marginal Likelihood.
% Each GP has a scale and bandwidth parameter that needs to be optimized. We
% could use different parameters for all or use the same parameter for all.
% This should be specified in hyperParams via fields useSameSm, useSamePr
% hyperParams should mandatorily contain the following fields
% - useSameSm, useSamePr: use the same smoothness/ scale values
% - fixSm, fixPr: fix the smoothness/ scale values - don't optimize over them.
% The following are optional
% - sigmaSmRange, sigmaPrRange: if using the same bw/scale the range for
%                               optimizing over them.
% - sigmaSmRanges, sigmaPrRanges: same as above, but if using different values
%     the ranges should be given here.
% This assumes that the decomposition is known and picks the other kernel hyper
% parameters.

  % prelims
  numDims = size(X, 1);
  numPts = size(X, 2);
  numGroups = numel(decomposition);
  useSameSm = hyperParams.useSameSm;
  useSamePr = hyperParams.useSamePr;
  diRectOptions.maxits = 8;

  % Set the hyperparameters for each GP
  % Common Mean Function 
  % -----------------------------------------
  if isempty(hyperParams.commonMeanFunc)
    commonMeanFunc = @(arg) mean(y) * ones(size(arg,1), 1);
%     commonMeanFunc = @(arg) zeros(size(arg,1), 1);
  else
    commonMeanFunc = hyperParams.commonMeanFunc;
%     commonMeanFunc,
  end
  % Common Noise parameter
  % -------------------------------------
    commonNoise = hyperParams.commonNoise;
  % Mean Functions for each GP
  % -------------------------------------
    if isempty(hyperParams.meanFuncs)
      hyperParams.meanFuncs = @(arg) zeros( size(arg,1), 1);
    end
    if numel(hyperParams.meanFuncs) == 1
      meanFuncs = cell(numGroups, 1);
      [meanFuncs{:}] = deal(hyperParams.meanFuncs);
    else
      meanFuncs = hyperParams.meanFuncs;
    end
  % Noise Parameters
  % -------------------------------------
    if numel(hyperParams.noises) == 1
      noises = hyperParams.noises * ones(numGroups, 1);
    else
      noises = hyperParams.noises;
    end

  % Prepare bounds for the parameters. 
  % -------------------------------------
  % If using the same bw across all groups
  if useSameSm
    if isempty(hyperParams.sigmaSmRange)
      sigmaSmRange = norm(std(X));
    else
      sigmaSmRange = hyperParams.sigmaSmRange;
    end
    % Set the Bounds for Optimization
    if size(sigmaSmRange, 2) == 2, sigmaSmBound = sigmaSmRange;
    else, sigmaSmBound = sigmaSmRange * [0.01 100];
    end
  else
  % Different Bandwidths for the different groups
    if isempty(hyperParams.sigmaSmRanges)
      sigmaSmRanges = zeros(numGroups, 1);
      for k = 1:numGroups
        sigmaSmRanges(k) = norm( std( X(:, decomposition{k}) ) );
      end
    else
      sigmaSmRanges = hypeParams.sigmaSmRanges;
    end
    % Set the bounds for optimization
    if size(sigmaSmRanges, 2) == 2, sigmaSmBounds = sigmaSmRanges;
    else sigmaSmBounds = sigmaSmRanges * [0.01 100];
    end
  end

  % Now do the same for the scale parameter
  if useSamePr
    if isempty(hyperParams.sigmaPrRange)
      sigmaPrRange = std(y);
    else
      sigmaPrRange = hyperParams.sigmaPrRange;
    end
    % Set the Bounds for Optimization
    if size(sigmaPrRange, 2) == 2, sigmaPrBound = sigmaPrRange;
    else, sigmaPrBound = sigmaPrRange * [0.1 10];
    end
  else
  % Different Scales for the different groups
    if isempty(hyperParams.sigmaPrRanges)
      sigmaPrRanges = std(y) * ones(numGroups, 1);
    else
      sigmaPrRanges = hyperParams.sigmaPrRanges;
    end
    % Set the Bounds for Optimization
    if size(sigmaPrRanges, 2) == 2, sigmaPrBounds = sigmaPrRanges;
    else, sigmaPrBounds = sigmaPrRanges * [0.1 10];
    end
  end

  % Not adding noise here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Since there are several noise parameters, its best if the caller explicitly
  % sets the noise parametrs.

  % If fixing the smoothness/ scale, set the optimal values right away
  if hyperParams.fixSm
    if useSameSm, sigmaSmOpts = sigmaSmRange * ones(numGroups, 1);
    else, sigmaSmOpts = sigmaSmRange
    end
  end
  if hyperParams.fixPr
    if useSamePr, sigmaPrOpts = sigmaPrRange * ones(numGroups, 1);
    else, sigmaPrOpts = sigmaPrRange;
    end
  end

  % Now we need to do the optimization over different scenrios
  % ============================================================================
  % Prelims
  oneVec = ones(numGroups, 1);

  % 1. If fixing both parameters, nothing to be done
  % ------------------------------------------------
  if hyperParams.fixSm && hyperParams.fixPr

  % 2. If fixing only the Smoothness
  % ------------------------------------------------
  elseif hyperParams.fixSm
    % only optimize over the scale paraemeter 
    if useSamePr
      nlmlF = @(t) normMargLikelihood(sigmaSmOpts, exp(t) * oneVec, ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      [~, logOptPr] = diRectWrap(nlmlF, log(sigmaPrBound), diRectOptions);
      sigmaPrOpts = exp(logOptPr) * oneVec;
    else
      nlmlF = @(t) normMargLikelihood(sigmaSmOpts, exp(t), ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      [~, logSigmaPrOpts] = diRectWrap(nlmlF, log(sigmaPrBounds),diRectOptions);
      sigmaPrOpts = exp(logSigmaPrOpts);
    end

  % 3. If fixing only the Scale 
  % ------------------------------------------------
  elseif hyperParams.fixPr
    % only optimize over the bandwidth parameter
    if useSamePr
      nlmlF = @(t) normMargLikelihood(exp(t) * oneVec, sigmaPrOpts,  ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      [~, logOptSm] = diRectWrap(nlmlF, log(sigmaSmBound), diRectOptions);
      sigmaSmOpts = exp(logOptSm) * oneVec;
    else
      nlmlF = @(t) normMargLikelihood(exp(t), sigmaPrOpts, ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      [~, logSigmaSmOpts] = diRectWrap(nlmlF, log(sigmaSmBounds), diRectOptions);
      sigmaSmOpts = exp(logSigmaSmOpts);
    end

  % 4. If fixing neither
  % ------------------------------------------------
  else
    % 4.1 If using the same for both scale and bandwidth for all groups
    if useSameSm && useSamePr
      nlmlF = @(t) normMargLikelihood( exp(t(1))*oneVec, exp(t(2))*oneVec, ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      diRectBounds = log([sigmaSmBound; sigmaPrBound]);
      [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
      sigmaSmOpts = exp(optParams(1)) * oneVec;
      sigmaPrOpts = exp(optParams(2)) * oneVec;

    % 4.2 If using the same smoothness for all groups
    elseif useSameSm
      nlmlF = @(t) normMargLikelihood( exp(t(1))*oneVec, exp(t(2:end)) * oneVec, ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      diRectBounds = log([sigmaSmBound; sigmaPrBounds]);
      [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
      sigmaSmOpts = exp(optParams(1)) * oneVec;
      sigmaPrOpts = exp(optParams(2:end));

    % 4.3 If using the same scale
    elseif useSamePr
      nlmlF = @(t) normMargLikelihood( exp(t(1:numGroups)), exp(t(end)) * oneVec, ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      diRectBounds = log([sigmaSmBounds; sigmaPrBound]);
      [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
      sigmaSmOpts = exp(optParams(1:numGroups)); 
      sigmaPrOpts = exp(optParams(end)) * oneVec;

    % 4.4 If using different scales and bws
    else
      nlmlF = @(t) normMargLikelihood( exp(t(1:numGroups)), exp(t(numGroups+1:end)), ...
        decomposition, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
      diRectBounds = log([sigmaSmBounds; sigmaPrBounds]);
      [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
      sigmaSmOpts = exp(optParams(1:numGroups));
      sigmaPrOpts = exp(optParams(numGroups+1 : end)); 

    end %  if useSameSm && useSamePr

  end % if fixSm & fixPr
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % phew! that was crazy ...

  % Finally train the GP
  gpHyperParams.commonMeanFunc = commonMeanFunc;
  gpHyperParams.meanFuncs = meanFuncs;
  gpHyperParams.commonNoise = commonNoise;
  gpHyperParams.noises = noises;
  gpHyperParams.sigmaSms = sigmaSmOpts;
  gpHyperParams.sigmaPrs = sigmaPrOpts;
  [mu, KPost, Mus, KPosts, combinedFuncH, funcHs] = ...
    addGPRegression(X, y, Xtest, decomposition, gpHyperParams);

end


function nlml = normMargLikelihood(sigmaSms, sigmaPrs, decomposition, ...
  X, y, meanFuncs, commonMeanFunc, noises, commonNoise)
% Computes the normalized log likelihood of the model.

  % prelims
  numData = size(X, 1);
  Ky = combinedKernelNoise(X, decomposition, sigmaSms, sigmaPrs, noises, ...
                           commonNoise);
  L = stableCholesky(Ky);
  y_ = y - combinedMeanFunc(X, commonMeanFunc, meanFuncs, decomposition);
  alpha = L' \ (L \ y_);
  nlml = -1/2 * y_' * alpha - sum(log(diag(L))) - numData/2 * log(2*pi);
end

