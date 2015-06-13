function [maxVal, maxPt, boQueries, boVals, history] = addGPBO(...
  oracle, decomp, bounds, numIters, params)
% Implements Gaussian Process Bandits/ Bayesian Optimisation using Additive 
% Gaussian Processes. Please see our ICML 2015 paper: "High Dimensional Bayesian
% Optimisation and Bandits via Additive Models". K.Kandasamy, J.Schneider, B.Poczos
%
% Inputs:
%   oracle: A function handle for the function you wish to optimise.
%   decomp: A cell/structure specifying the decomposition (see below).
%   bounds: A numDimsx2 matrix specifying the lower and upper bound of each
%   dimension.
%   numIters: The number of GPB/BO iterations
%   params: A structure specifying the various hyper parameters for optimisation. If
%     you wish to use the default settings, pass an empty struct. Otherwise, see
%     demoCustomise.m to see how to set each hyper parameter. Also see below.
%
% Outputs:
%   maxVal: The maximum queried value of the function.
%   maxPt: The queried point with the maximum value.
%   boQueries: A matrix indicating the points at which the algorithm queried.
%   boVals: A vector of the query values.
%   history: A vector of the maximum value obtained up until that iteration.
%
% params should have a field called decompStrategy: It should be one of 
% 'known', 'learn', 'random' and 'partialLearn'. 
% 'known': The decomposition is known and given in decomp. We will optimise
% according to this.
% For the other 3 cases, decomp should have two parameters d & M.
% 'learn': The decomposition is unknown and should be learned.
% 'random': Randomly pick a partition at each iteration. 
% 'partialLearn': Partially learn the decomposition at each iteration by trying
%   out a few and picking the best.
% The default is partialLearn and is the best option if you don't know much about
% your function.

  % Define these to avoid typos
  DECOMP_KNOWN = 'known';
  DECOMP_LEARN = 'learn';
  DECOMP_RAND = 'random';
  DECOMP_PLEARN = 'partialLearn';

  % Prelims
  numDims = size(bounds, 1);
  dummyPts = zeros(0, numDims); % to build the GPs
  MAX_THRESHOLD_EXCEEDS = 5;
  NUM_ITERS_PER_PARAM_RELEARN = 25;

  % First check for all parameters and if not given use default values
  if ~exist('params', 'var') 
    params = struct();
  end
  % Determine the number of groups
  if isfield(decomp, 'M'), numGroups = decomp.M;
  else numGroups = numel(decomp);
  end
  % Screen for hyper parameters and otherwise set default values.
  params = boProcessParams(params, oracle, bounds, numGroups);

  % The Decomposition
  % -----------------------------------------------------
  gpHyperParams.decompStrategy = params.decompStrategy;
    if strcmp(params.decompStrategy, DECOMP_KNOWN) 
      decomposition = decomp;
      numGroups = numel(decomposition);
      % Do some diagnostics on the decomposition and print them out
      relevantCoords = cell2mat(decomposition');
      numRelevantCoords = numel(relevantCoords);
      if numRelevantCoords ~= numel(unique(relevantCoords))
        error('The same coordinate cannot appear in different groups');
      end
      fprintf('# Groups: %d, %d/%d coordinates are relevant\n', ...
        numGroups, numRelevantCoords, numDims);

    elseif isfield(decomp, 'M')
      % Now decomposition should have two fields d and M
      numGroups = decomp.M;

    else % in this case the decomposition is given.
      numGroups = numel(decomp);
    end

  % Initialisation points
  initPts = params.initPts;
  initVals = params.initVals;
  numInitPts = size(initPts, 1);

  % Hyper parameters for the GP
  gpHyperParams.decompStrategy = params.decompStrategy;
  gpHyperParams.commonMeanFunc = params.commonMeanFunc;
  gpHyperParams.meanFuncs = params.meanFuncs;
  gpHyperParams.commonNoise = params.commonNoise;
  gpHyperParams.noises = params.noises;
  gpHyperParams.fixPr = params.fixPr;
  gpHyperParams.useSamePr = params.useSamePr;
  gpHyperParams.sigmaPrRange = params.sigmaPrRange;
  gpHyperParams.sigmaPrRanges = params.sigmaPrRanges;
  gpHyperParams.useSameSm = params.useSameSm;

  % The Bandwidth 
  if params.useFixedBandWidth
    gpHyperParams.fixSm = true;
    if params.useSameSm, gpHyperParams.sigmaSmRange = params.sigmaSmRange;
    else, gpHyperParams.sigmaSmRanges = params.sigmaSmRanges;
    end
  else % This BO algorithm will set the bw via its own procedure
    alBWLB = params.alBWLB;
    alBWUB = params.alBWUB;
    % Set an initial bw. This will change as the algorithm progresses
    alCurrBW = alBWUB;
  end

  % Define the following before proceeding
  boQueries = [initPts; zeros(numIters, numDims)];
  boVals = [initVals; zeros(numIters, 1)];
  history = [max(initVals(cumsum(triu(ones(length(initVals))))))'; ...
             zeros(numIters, 1)];
  threshExceededCounter = 0;
  if isempty(initVals)
    currMaxVal = -inf;
    currMaxPt = [];
  else
    [currMaxVal, maxIdx] = max(initVals);
    currMaxPt = initPts(maxIdx, :);
  end
  
  fprintf('Peforming BO (dim = %d)\n', numDims);
  for boIter = 1:numIters

    if mod(boIter, NUM_ITERS_PER_PARAM_RELEARN) == 0
      fprintf('Additive GP BO iter %d/ %d. MaxVal: %0.4f CumReward: %0.4f\n',...
        boIter, numIters, currMaxVal, sum(boVals)/(boIter+numInitPts) );
    end

    % Prelims
    numBoPts = numInitPts + boIter - 1;
    currBoQueries = boQueries(1:numBoPts, :);
    currBoVals = boVals(1:numBoPts);

    % First redefine ranges for the GP bandwidth if needed
    if ~params.useFixedBandWidth & ...
       ( mod(boIter-1, NUM_ITERS_PER_PARAM_RELEARN) == 0 | ...
         threshExceededCounter == MAX_THRESHOLD_EXCEEDS )

      if threshExceededCounter == MAX_THRESHOLD_EXCEEDS
        alBWUB = max(alBWLB, 0.9 * alCurrBW);
        threshExceededCounter = 0;
        fprintf('Threshold Exceeded %d times - Reducing BW\n', ...
                MAX_THRESHOLD_EXCEEDS);
      else
%         alBWUB = max(alBWLB, 0.9 * alBWUB);
      end

      % Define the BW range for addGPMargLikelihood      
      if alBWUB == alBWLB
        gpHyperParams.fixSm = true;
        gpHyperParams.sigmaSmRanges = alBWLB * ones(numGroups, 1);

      else
        gpHyperParams.fixSm = false;
        % Use same bandwidth for now. TODO: modify to allow different bandwidths
        gpHyperParams.useSameSm = true;
        gpHyperParams.sigmaSmRange = [alBWLB, alBWUB];
      end % alBWUB == alBWLB

      % Obtain the optimal GP parameters
      [~, ~, ~, ~, ~, ~, ~, alCurrBWs, alCurrScales, ~, learnedDecomp] = ...
        addGPDecompMargLikelihood( currBoQueries, currBoVals, dummyPts, ...
          decomp, gpHyperParams );

      alCurrBW = alCurrBWs(1); %TODO: modify to allow different bandwidths
      if numDims < 24
        learnedDecompStr = '';
        for i = 1:numel(learnedDecomp)
          learnedDecompStr = [learnedDecompStr mat2str(learnedDecomp{i})];
        end
      else
        learnedDecompStr = '';
      end
      fprintf('Picked bw: %0.4f (%0.4f, %0.4f), Scale: %0.4f. Decomp: %s (%d)\n', ...
        alCurrBW, alBWLB, alBWUB, alCurrScales(1), ...
        learnedDecompStr, numel(learnedDecomp) );

    end % if ~params.useFixedBandWidth ...

    % Now build the GP
    currGPParams.commonMeanFunc = gpHyperParams.commonMeanFunc;
    currGPParams.meanFuncs = gpHyperParams.meanFuncs;
    currGPParams.commonNoise = gpHyperParams.commonNoise;
    currGPParams.noises = gpHyperParams.noises;
    currGPParams.sigmaSms = alCurrBWs;
    currGPParams.sigmaPrs = alCurrScales;
    [~,~,~,~, combFuncH, funcHs] = addGPRegression(currBoQueries, ...
      currBoVals, dummyPts, learnedDecomp, currGPParams);

    % Now obtain the next point
    [nextPt, ~, nextPtStd] = getNextQueryPt(params, combFuncH, funcHs, ...
      learnedDecomp, currBoVals, bounds);
    % If it is too close, perturb it a bit
    if min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2, 2) ))/...
          alCurrBW< 1e-10 
      while min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2,2)))/...
                 alCurrBW < 1e-10
        nextPt = projectToRectangle( ...
          nextPt' + 0.1 * alCurrBW * randn(numDims, 1), bounds)';
      end
    end

    % Determine the current best point
    nextPtVal = oracle(nextPt);
    if nextPtVal > currMaxVal
      currMaxVal = nextPtVal;
      currMaxPt = nextPt;
    end
    boQueries(numInitPts + boIter, :) = nextPt;
    boVals(numInitPts + boIter) = nextPtVal;
%     fprintf('#: %d, maxVal: %0.5f, currVal: %0.5f\n', ...
%       boIter, currMaxVal, nextPtVal);

    % Check if nextPtStd is too small
    if nextPtStd < params.optPtStdThreshold
      threshExceededCounter = threshExceededCounter + 1;
    else
      threshExceededCounter = 0;
    end

    % Store the current best value
    history(boIter+numInitPts) = currMaxVal;

  end % end for boIter

  maxVal = currMaxVal;
  maxPt = currMaxPt;

end


% function to obtain the next query point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nextPt, nextPtMean, nextPtStd, nextPtUtil] = ...
  getNextQueryPt(params, combFuncH, funcHs, decomposition, boVals, bounds)

  % So here's the plan here. We will define a utility function over each GP and
  % pick the point that maximizes this utility function
  
  numGroups = numel(decomposition);
  numDims = size(bounds, 1);

  nextPt = zeros(1, numDims);
  nextPtUtils = zeros(numGroups, 1);

  for k = 1:numGroups
    
    coords = decomposition{k};
    currBounds = bounds(coords, :);
    currGPFuncH = funcHs{k};

    if strcmp(params.utilityFunc, 'UCB')
      utility = @(t) getUCBUtility(t, currGPFuncH, size(boVals, 1)); 
    elseif strcmp(params.utilityFunc, 'EI')
      % Only applies to non-additive models. (i.e. d = D, M = 1)
      utility = @(t) getEIUtility(t, currGPFuncH, max(boVals));
    else
      utility = @(t) params.utilityFunc(t, currGPFuncH);
    end

    % Run DiRect
    [currNextPtUtil, currNextPt, hist] = ...
      diRectWrap(utility, currBounds, params.diRectParams);

    % store nextPt in relevant coordinates
    nextPt(coords) = currNextPt;
    nextPtUtils(k) = currNextPtUtil;

  end

  % Finally return the mean and the standard deviation at nextPt
  [nextPtMean, nextPtStd] = combFuncH(nextPt);

end

