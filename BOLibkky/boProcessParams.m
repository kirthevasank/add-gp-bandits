function params = boProcessParams(params, oracle, bounds, numGroups)
% Check for parameters and if not use default values.

  % Prelims
  numDims = size(bounds, 1);

  % Initialisation
  if ~isfield(params, 'initPts') | isempty(params.initPts)
    if ~isfield(params, 'numInitPts')
      params.numInitPts = 10;
    end
    fprintf('Obtaining %d Points for Initialisation.\n', params.numInitPts);
    params.initPts = boGetInitPts(bounds, params.numInitPts);
    params.initVals = oracle(params.initPts);
  else
    params.numInitPts = size(params.initPts, 1);
    fprintf('Using the given %d Initialisation points.\n', params.numInitPts);
  end
  initPts = params.initPts;
  initVals = params.initVals;
  stdInitVals = std(initVals);

  % Hyper Parameters for optimisation
  % ============================================================================
  % Utility (Acquisition) Function - use UCB as default
  if ~isfield(params, 'utilityFunc') | isempty(params.utilityFunc)
    params.utilityFunc = 'UCB';
  end
  % Number of DiRect evaluations to optimise the utility function
  if ~isfield(params, 'numDiRectEvals') | isempty(params.numDiRectEvals)
    params.numDiRectEvals = min(5000, 100*numDims);
  end
  % The not-so-important hyper parameters
  if ~isfield(params, 'optPtStdThreshold') | isempty(params.optPtStdThreshold)
    params.optPtStdThreshold = 0.002;
  end

  % Set up sub-structure for DiRect
  diRectParams.maxevals = ceil(params.numDiRectEvals/numGroups);
  diRectParams.maxits = inf;
  params.diRectParams = diRectParams;

  % GP Hyper parameters for Regression
  % ============================================================================
  % The Common Mean Function
  if ~isfield(params, 'commonMeanFunc')
    if isfield(params, 'meanFuncValue')
      params.commonMeanFunc = ...
        @(arg) params.meanFuncValue * ones(size(arg,1), 1);
    else
      params.commonMeanFunc = []; % By default will use all zeros.
    end
  end
  % The Mean function for the individual GPs
  if ~isfield(params, 'meanFuncs') | isempty(params.meanFuncs)
    params.meanFuncs = @(arg) zeros( size(arg,1), 1);
  end
  % Common noise parameters
  if ~isfield(params, 'commonNoise') | isempty(params.commonNoise)
    params.commonNoise = 0.01 * stdInitVals;
  end
  % Individual noise
  if ~isfield(params, 'noises') | isempty(params.noises)
    params.noises = 0;
  end
  % Scale parameters
  % ----------------
  if ~isfield(params, 'fixPr')
    params.fixPr = false;
  end
  if ~isfield(params, 'useSamePr')
    params.useSamePr = true;
  end
  if ~isfield(params, 'sigmaPrRange')
    params.sigmaPrRange = [0.03 30] * stdInitVals;
  end
  if ~isfield(params, 'sigmaPrRanges')
    params.sigmaPrRanges = []; % use same prior by default
  end
  % Bandwidth parameters
  % --------------------
  if ~isfield(params, 'useFixedBandWidth')
    params.useFixedBandWidth = false;
  end
  if ~isfield(params, 'fixSm')
    params.fixSm = false;
  end
  if ~isfield(params, 'useSameSm')
    params.useSameSm = true;
  end
  if ~isfield(params, 'alBWLB')
    params.alBWLB = 1e-5;
  end
  if ~isfield(params, 'alBWUB')
    params.alBWUB = 5;
  end
    


end
