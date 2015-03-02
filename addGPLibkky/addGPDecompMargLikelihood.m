function [mu, KPost, Mus, KPosts, combinedXFuncH, combinedZFuncH, funcHs, ...
  sigmaSmOpts, sigmaPrOpts, A, learnedDecomp] = ...
  addGPRotDecompLikelihood(X, y, Xtest, decomp, hyperParams)
% This function attempts to find the best Kernel hyper parameters to fit an
% additive function. The kernel parameters include the smoothness and scale
% parameters and the decomposition.
% Inputs:
%   X, y, Xtest: Training and Testing data
%   decomp: If the decomposition need not be learned, then this should contain
%     the true decomposition. Otherwise, it should contain two fields d (# of
%     dimensions per group) and M (number of groups). 
%   hyperParams should have a field called decompStrategy: It should be one of 
%   'known', 'learn', 'random' and 'partialLearn'. 
%   'known': The decomposition is known and given in decomp. We will optimize
%   according to this.
%   For the other 3 cases, decomp should have two parameters d & M.
%   'learn': The decomposition is unknown and should be learned.
%   'random': Randomly pick a partition at each iteration. 
%   'partialLearn': Partially learn the decomposition at each iteration by trying
%     out a few and picking the best.

  % Define these to avoid typos
  DECOMP_KNOWN = 'known';
  DECOMP_LEARN = 'learn';
  DECOMP_RAND = 'random';
  DECOMP_PLEARN = 'partialLearn';

  % prelims
  D = size(X, 2);
  numDims = D;
  numPts = size(X, 1);
  % for diRect
  diRectOptions.maxits = 8;

  if ~strcmp(hyperParams.decompStrategy, DECOMP_KNOWN)

    if isfield(decomp, 'M')
      % First create a placeholder for the decomposition
      M = decomp.M;
      d = decomp.d;
      p = d*M;
      decomposition = cell(M, 1);
      for i = 1:M
        decomposition{i} = ((i-1)*d +1) : (i*d) ;
      end
    else % then decomp is a vector of values with the number of dims in each group.
      M = numel(decomp);
      d = max(decomp);
      p = sum(decomp);
      cumDims = cumsum(decomp); cumDims = cumDims(:);
      cumDims = [0; cumDims];
      decomposition = cell(M, 1);
      for i = 1:M
        decomposition{i} = [(cumDims(i)+1):cumDims(i+1)];
      end
    end
  else
    decomposition = decomp;
  end
  numGroups = numel(decomposition);
  oneVec = ones(numGroups, 1);

  if ~isfield(hyperParams, 'numBwSigmaDiRectIters') | ...
    isempty(hyperParams.numBwSigmaDiRectIters)
    numBwSigmaDiRectIters = 20;
  else numBwSigmaDiRectIters = hyperParams.numBwSigmaDiRectIters;
  end

  % Set the Hyperparameters for each GP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

  % Some parameters for learning the Kernel smoothness and range
  % Define Bounds for optimization of the bandwidth and scale
  % ---------------------------------------------------------
  sigmaSmBound = hyperParams.sigmaSmRange;
  sigmaPrBound = hyperParams.sigmaPrRange;
  diRectOptions.maxits = numBwSigmaDiRectIters;
  diRectBounds = log([sigmaSmBound; sigmaPrBound]);

  % Learn the Decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~strcmp(hyperParams.decompStrategy, DECOMP_KNOWN)

    % As the first A, use an arbitrary orthogonal matrix
    A = randn(D, p); A = orth(A);

    nlmlF = @(t) normRotMargLikelihood( exp(t(1))*oneVec, exp(t(2))*oneVec, ...
      decomposition, A, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
    [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
    sigmaSmOpts = exp(optParams(1)) * oneVec;
    sigmaPrOpts = exp(optParams(2)) * oneVec;

    % Now optimze w.r.t A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    negMLF = @(T) negNormRotMargLikelihood( sigmaSmOpts, sigmaPrOpts, ...
      decomposition, T, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
    if strcmp(hyperParams.decompStrategy, DECOMP_RAND)
      A = getRandPermMat(D);
    elseif strcmp(hyperParams.decompStrategy, DECOMP_PLEARN)
      A = decompOptPartial(negMLF, D, d, M);
    elseif strcmp(hyperParams.decompStrategy, DECOMP_LEARN)
      A = decompOptBrute(negMLF, D, d, M);
    else
      error('Unknown Strategy to handle decomposition\n');
    end 

  else
  % Otherwise, use the given decomposition
    A = eye(D);
  end
  % Learn the Decomposition Ends here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Finally optimize w.r.t h and sigma again
  nlmlF = @(t) normRotMargLikelihood( exp(t(1))*oneVec, exp(t(2))*oneVec, ...
    decomposition, A, X, y, meanFuncs, commonMeanFunc, noises, commonNoise);
  [~, optParams] = diRectWrap(nlmlF, diRectBounds, diRectOptions);
  sigmaSmOpts = exp(optParams(1)) * oneVec;
  sigmaPrOpts = exp(optParams(2)) * oneVec;

  % Finally Train the GP
  gpHyperParams.commonMeanFunc = commonMeanFunc;
  gpHyperParams.meanFuncs = meanFuncs;
  gpHyperParams.commonNoise = commonNoise;
  gpHyperParams.noises = noises;
  gpHyperParams.sigmaSms = sigmaSmOpts;
  gpHyperParams.sigmaPrs = sigmaPrOpts;
  % The function handles returned are for the transformed space but that is what
  % we want.
  Z = X * A;
  [mu, KPost, Mus, KPosts, combinedZFuncH, funcHs] = ...
    addGPRegression(Z, y, Xtest, decomposition, gpHyperParams);
  combinedXFuncH = @(X) combinedZFuncH(X*A);

  % Finally return the learned decomposition
  [~, permutedOrder] = orthToPermutation(A);
  learnedDecomp = cell(numGroups, 1);
  for i = 1:numGroups
    learnedDecomp{i} = permutedOrder(decomposition{i});
  end

end

