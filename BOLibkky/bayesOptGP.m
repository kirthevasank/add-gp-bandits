function [maxVal, maxPt, boQueries, boVals, history] = bayesOptGP(oracle, ...
  bounds, numIters, params)
% Performs Bayesian Optimization to maximize the function specified in the
% handle 'func'.
% initPts and initVals are points at which the function has been
% evaluated already.
% bounds are the bounds of the parameter space.
% numIters is the number of Iterations
% params contains the following fields % - utilityFunc: Give your own criterion
% here or use an already available
%    utility funciton (so far just 'EI' for Expected Improvement
% - optPtStdThreshold: The t_sigma param in (Wang, NdeF et al, IJCAI'13). If the
%   variance of the chosen point is consecutively smaller you should decrease
%   the GP bandwidth.
% - alBWLB, alBWUB: Lower and Upper Bounds for the BO GP bandwidth. Applies only
%   if alBandWidth (see below) is not specified.
% - numInitPts: number of points to initialize the problem.
% - alBandWidth: the bandwidth for the BO GP
% - gpNoiseLevel: the Noise Level for GPs
% - logLiklRange: the range for the likelihood. Something in the ballpark is ok.
% - meanFuncValue: the value for the mean function of the GP. Sometimes setting
%                  this to be a low value may be useful. 
% The other parameters of the GP are picked via marginal likelihood criterion
% If params.useFixedBandwidth is false, then params.postVarThreshold should be
% specified.

  % Prelims
  numDims = size(bounds, 1);
  dummyPts = zeros(0, numDims); % to build GPs
  MAX_THRESHOLD_EXCEEDS = 5; %# times the variance is allowed to fall low
  NUM_ITERS_PER_PARAM_RELEARN = 20;

  % If Init points are not given, initialize
  if ~isfield(params, 'initPts') | isempty(params.initPts)
    initPts = boGetInitPts(bounds, params.numInitPts);
    initVals = oracle(initPts);
  else
    initPts = params.initPts;
    initVals = params.initVals;
  end
  numInitPts = size(initPts, 1);

  % Check for parameters expected in params
  if ~isfield(params, 'diRectParams')
    params.diRectParams = struct;
    fprintf('WARNING: diRectParams not given. DiRect will use default vals.\n');
    params.diRectParams.maxits = 100;
    params.diRectParams.maevals = min(10000, 2^numDims);
  end
  if ~isfield(params, 'useFixedBandwidth')
    params.useFixedBandwidth = false;
  end

  % Determine values to be sent to the GP
  % =========================================================
  % The Mean Function
  % --------------------------------------
    if isfield(params, 'meanFunc')
      gpHyperParams.meanFunc = params.meanFunc; 
    else 
      if isfield(params, 'meanFuncValue')
        gpHyperParams.meanFunc = ...
          @(arg) params.meanFuncValue * ones(size(arg,1), 1);
      else
        gpHyperParams.meanFunc = []; 
      end 
    end 
  % The Scale
  % --------------------------------------
    if isfield(params, 'sigmaPr')
      % If the parameters specify the scale to use, then don't look beyond
      alCurrScale = params.sigmaPr;
    elseif isfield(params, 'scaleRange')
      alCurrScale = std(initVals);
      gpHyperParams.sigmaPrRange = params.scaleRange;
      gpHyperParams.sigmaPr =  0; % if we set this
    else ~isfield(params, 'scaleRange')
      % If the scale isn't given, but a range is given then we need to do
      % MargLikelihood. Specify this in the params.
      % If not given, then GPLib will use the regressands to estimate this.
      alCurrScale = std(initVals);
      gpHyperParams.sigmaPrRange = []; 
      gpHyperParams.sigmaPr =  0; % if we set this
    end
  % The Bandwidth
  % --------------------------------------
    if params.useFixedBandwidth
    % If use a fixed bandwidth, set it now. Otherwise we will set at each iter.
      gpHyperParams.sigmaSm = params.alBandWidth;
    else
      alBWLB = params.alBWLB;
      alBWUB = params.alBWUB;
      % Set an initial BW. This will change as the algorithm progresses.
      alCurrBW = alBWUB; 
    end
  % The Noise Level 
  % --------------------------------------
  if ~isfield(params, 'gpNoiseLevel')
    gpHyperParams.noise = 0.01 * std(initVals);
  else
    gpHyperParams.noise = params.gpNoiseLevel;
  end

  % Define the following before proceeding
  boQueries = [initPts; zeros(numIters, numDims)];
  boVals = [initVals; zeros(numIters, 1)];
  history = [max(initVals(cumsum(triu(ones(length(initVals))))))'; ...
             zeros(numIters, 1)];
  threshExceededCounter = 0;
  [currMaxVal, maxIdx] = max(initVals);
  currMaxPt = initPts(maxIdx, :);
  
  fprintf('Performing BO (dim = %d)\n', numDims);
  for boIter = 1:numIters

    if mod(boIter, 20) == 0
      fprintf('Bayesian Optimization iter %d/ %d. MaxVal: %0.4f\n', ...
        boIter, numIters, currMaxVal);
    end
    % prelims
    numBoPts = numInitPts + boIter - 1;
    currBoQueries = boQueries(1:numBoPts, :);
    currBoVals = boVals(1:numBoPts);

    % First rebuild the GP if needed
    if ~params.useFixedBandwidth
      if mod(boIter-1, NUM_ITERS_PER_PARAM_RELEARN) == 0 | ...
          threshExceededCounter == MAX_THRESHOLD_EXCEEDS
        if threshExceededCounter == MAX_THRESHOLD_EXCEEDS
          alBWUB = max(alBWLB, 0.9 * alCurrBW);
          threshExceededCounter = 0;
          fprintf('Threshold exceeded %d times- Reducing Bandwidth\n', ...
                   MAX_THRESHOLD_EXCEEDS);
        else
          alBWUB = max(alBWLB, 0.9 * alBWUB);
        end

        % Learn the optimal parameters for the GP Model.
        if alBWUB == alBWLB
          gpHyperParams.sigmaSm = alBWLB;
        else
          gpHyperParams.sigmaSm = 0;
          gpHyperParams.sigmaSmRange = [alBWLB, alBWUB];
        end
        [~, ~, ~, alCurrBW, alCurrScale] = GPMargLikelihood(currBoQueries, ...
          currBoVals, dummyPts, gpHyperParams);        
        % See end of page for some debugging code here.
      end
    end

    % Build the GP: returns the function handle
    runTimeParams.retFunc = true;
    currGPParams.meanFunc = gpHyperParams.meanFunc;
    currGPParams.sigmaSm = alCurrBW;
    currGPParams.sigmaPr = alCurrScale;
    currGPParams.noise = gpHyperParams.noise;
    [~, ~, ~, gpFuncH] = GPRegression(currBoQueries, currBoVals, dummyPts, ...
                           currGPParams, runTimeParams);

    % Obtain the next query point and query
    [nextPt, ~, nextPtStd] = getNextQueryPt(params, gpFuncH, currBoVals, bounds);
    % If it is too close, then perturb it a bit
    if min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2, 2) ))/alCurrBW < 1e-10 
%       threshExceededCounter = threshExceededCounter + 1;
      while min( sqrt( sum( bsxfun(@minus, currBoQueries, nextPt).^2, 2)))/alCurrBW < 1e-10
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

  end % for boIter

  maxVal = currMaxVal;
  maxPt = currMaxPt;
%   [maxVal, maxIdx] = max(boVals);
%   maxPt = boQueries(maxIdx, :);
    
end


function [nextPt, nextPtMean, nextPtStd, nextPtUtil] = ...
  getNextQueryPt(params, gpFuncH, boVals, bounds)
% This is what this function should do. It should maximize the utility 

  % First create the utility function to be maximized
  if strcmp(params.utilityFunc, 'EI')
    utility = @(t) getEIUtility(t, gpFuncH, max(boVals));
  elseif strcmp(params.utilityFunc, 'UCB')
    utility = @(t) getUCBUtility(t, gpFuncH, size(boVals, 1));
  else
    utility = @(t) params.utilityFunc(t, gpFuncH);
  end

  % Now run DiRect
  [nextPtUtil, nextPt, hist] = diRectWrap(utility, bounds, params.diRectParams);

  % Now get the mean and variance of nextPt
  [nextPtMean, nextPtStd] = gpFuncH(nextPt);

end



% DEBUG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num_samples = 100;
% N = 100; Z = linspace(0,1,N)';
% [post_mean, K, funcH, sigmaSm, sigmaPr] = ...
%   GPMargLikelihood(currBoQueries, currBoVals, Z, gpHyperParams);
% 
% %   % Now draw some samples
% %   num_samples = 100;
% %   gp_samples = GPDrawSamples(post_mean, K, num_samples);
% %   % plot the samples
% %   figure;
% %   hold on,
% %   for i = 1:num_samples
% %     plot(Z, gp_samples(i,:), 'm-');
% %   end
% plot(Z, oracle(Z), 'g', 'LineWidth', 4); hold on,
% plot(Z, post_mean, 'b', 'LineWidth', 2); hold on,
% plot(currBoQueries, currBoVals, 'kx', 'MarkerSize', 10, 'Linewidth', 2);
% % [currBoQueries, currBoVals],
% fprintf('Chosen: sigmaSm: %f, sigmaPr: %f, \nlikl of post-mean: %f\n', ...
%   sigmaSm, sigmaPr, GPAvgLogLikelihood(post_mean, K, post_mean));
% fprintf('Paused .... \n');
% pause;
% close;
