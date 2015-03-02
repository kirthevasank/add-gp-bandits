function [mu, KPost, Mus, KPosts, combinedFuncH, funcHs] = ...
  additiveGPRegression(X, y, Xtest, decomposition, hyperParams)
% A matlab function for perform GP Regression when the model is additive. This
% performs inference in each individual GP.
% Inputs
%   X, y, Xtest: Training input/output and test input. If no test input set
%                Xtest to a 0xnumDims vector
%   decomposition: A matlab cell array containing numGroups elements. Each
%                  element is a vector containing the coordinates in that group.
%   hyperParams: Contains the smoothness (sigmaSms), scale (sigmaPrs) and noise
%                (noise0, noises) parameters for each GP in the additive model. 
% Outputs
% Mus, KPosts: The predictive mean and covariance for each GP in the additive
%              model.
% addFuncH: A function handle for the combined GP.
% funcHs: A cell array with a function handle for each group.

  % Prelims
  numTrData = size(X, 1);
  numDims = size(X, 2);
  numGroups = numel(decomposition);

  % Set the hyperparameters for each GP
  % Common Mean Function 
  % -----------------------------------------
  if isempty(hyperParams.commonMeanFunc)
    commonMeanFunc = @(arg) mean(y) * ones(size(arg,1), 1);
  else
    commonMeanFunc = hyperParams.commonMeanFunc;
  end
  % Common Noise parameter
  % -------------------------------------
    commonNoise = hyperParams.commonNoise;
  % Smoothness parameters
  % -------------------------------------
    if numel(hyperParams.sigmaSms) == 1
      sigmaSms = hyperParams.sigmaSms * ones(numGroups, 1);
    else
      sigmaSms = hyperParams.sigmaSms;
    end
  % Scale parameters
  % -------------------------------------
    if numel(hyperParams.sigmaPrs) == 1
      sigmaPrs = hyperParams.sigmaPrs * ones(numGroups, 1);
    else
      sigmaPrs = hyperParams.sigmaPrs;
    end
  % Mean Functions for each GP
  % -------------------------------------
    if isempty(hyperParams.meanFuncs)
      hyperParams.meanFuncs = @(arg) zeros( size(arg,1), 1);
    end
    if numel(hyperParams.meanFuncs) == 1 && iscell(hyperParams.meanFuncs)
      hyperParams.meanFuncs = hyperParams.meanFuncs{1};
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

  if numTrData < inf
    
    % Construct the Training Kernel Matrix and Invert it
    % --------------------------------------------------
    K0 = combinedKernelNoise(X, decomposition, sigmaSms, sigmaPrs, ...
                             noises, commonNoise);
    % To Invert this we need to do the cholesky decomposition
    L = stableCholesky(K0); 
    % Compute alpha
    alpha = L'\ (L \ (y - ...
      combinedMeanFunc(X, commonMeanFunc, meanFuncs, decomposition)));

    % Now obtain function Handles for the individual GPs
    % --------------------------------------------------
    % First for the individual GP Outputs
    funcHs = cell(numGroups, 1);
    for k = 1:numGroups
      coords = decomposition{k};
      bw = sigmaSms(k);
      scale = sigmaPrs(k);
      currMeanFunc = meanFuncs{k};
      funcHs{k} = @(Xte) computeIndGPOutput(Xte, X, coords, L, alpha, ...
        bw, scale, currMeanFunc);
    end 
    % Now the combined GP Output
    combinedFuncH = @(Xte) computeCombGPOutput(Xte, X, decomposition, L, ...
      alpha, sigmaSms, sigmaPrs, meanFuncs, commonMeanFunc);

    % Now obtain the Ouputs
    % --------------------------------------------------
    % First the individual GP Outputs
      Mus = cell(numGroups, 1);
      KPosts = cell(numGroups, 1);
      for k = 1:numGroups
        [currMu, ~, currK] = funcHs{k}(Xtest);
        Mus{k} = currMu;
        KPosts{k} = currK;
      end
    % Now the combined GP output
    [mu, ~, KPost] = combinedFuncH(Xtest);

  else % numTrainData <
    % TODO: work on a more efficient implementation for large data.
  end % numTrainData <

end % end Main Function

