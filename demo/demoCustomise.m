% This shows how to tune the various hyper parameters in our implementation. We have
% used good heuristics as default so unless you really understand the code well we
% would recommend sticking to the defaults.

close all;
clear all;
addpath ../addGPLibkky/
addpath ../BOLibkky/
addpath ../utils/
warning off;


% Problem parameters
numExperiments = 2; % Running 2 experiments to generate the error bars.
numDims = 10; trueNumDimsPerGroup = 4; numDimsPerGroupCands = [10 1 2 5]';
% numDims = 24; trueNumDimsPerGroup = 12; numDimsPerGroupCands = [24 1 3 6 12]'; 
% numDims = 40; trueNumDimsPerGroup = 18; numDimsPerGroupCands = [40 1 5 10 20]';
% We will run multiple instantiations of add-gp-ucb with maximum group size given in
% numDimsPerGroupCands. In the first case (d = numDims) it is plain gp-ucb.

% Experiment parameters
numIters = 100;
numDiRectEvals = min(5000, 100*numDims);

numdCands = numel(numDimsPerGroupCands);
% Get the function
[func, funcProperties] = getAdditiveFunction(numDims, trueNumDimsPerGroup);
bounds = funcProperties.bounds;
trueDecomp = funcProperties.decomposition;
trueMaxVal = funcProperties.maxVal;
trueMaxPt = funcProperties.maxPt;
trueNumGroups = numel(trueDecomp);

% Compute some statistics to help with the initialization
% This is something you won't be able to do with an expensive function. In the
% default settings we set this using the initialisation points.
th = rand(1000, numDims); fth = func(th);
meanFth = mean(fth);
stdFth = std(fth);

% Parameters for additive Bayesian optimization
boParams.optPtStdThreshold = 0.002;
boParams.useFixedBandwidth = false; % Don't use same kernel at all iterations
boParams.numInitPts = 10; %20; % min(20, numDims);
boParams.commonNoise = 0.01 * stdFth;
boParams.utilityFunc = 'UCB';
boParams.meanFuncs = []; % No mean funcs for the individual GPs
boParams.commonMeanFunc = @(arg) zeros(size(arg, 1), 1);
boParams.useSamePr = true; % Use same bandwidth and scale for all GPs
boParams.useSameSm = true;
boParams.fixPr = false;
boParams.fixSm = false;
% The upper/ lower bounds for the length scale parameter (called sigmaPr in this
% implementation).
boParams.sigmaPrRange = [0.03 30] * stdFth;
% The upper/ lower bounds for the bandwidth for the GP. To be selected by
% maximising the GP marginal likelihood.
boParams.alBWLB = 1e-5; 
boParams.alBWUB = 5;

% From here on customize each parameters separately.
% Known Decomposition
boKDParams = boParams;
boKDParams.decompStrategy = 'known';
boKDParams.diRectParams.maxevals = ceil(numDiRectEvals/trueNumGroups);
boKDParams.diRectParams.maxits = inf;
% The rest - arbitrary decompositions
boAddParams = boParams;
boAddParams.decompStrategy = 'partialLearn';
boAddParams.diRectParams.maxits = inf;
% EI
boEIParams = boParams;
boEIParams.utilityFunc = 'EI';
boEIParams.diRectParams.maxevals = numDiRectEvals;
boEIParams.diRectParams.maxits = inf;
boEIParams.decompStrategy = 'known';
boEIDecomp = {[1:numDims]};

totalNumQueries = numIters + boParams.numInitPts;
% Initialize arrays for storing the history
boKDHistories = zeros(numExperiments, totalNumQueries);
boAddHistories = zeros(numExperiments, totalNumQueries, numdCands);
boEIHistories = zeros(numExperiments, totalNumQueries);
randHistories = zeros(numExperiments, totalNumQueries);
% For storing simple regret values
boKDSimpleRegrets = zeros(numExperiments, totalNumQueries);
boAddSimpleRegrets = zeros(numExperiments, totalNumQueries, numdCands);
boEISimpleRegrets = zeros(numExperiments, totalNumQueries);
randSimpleRegrets = zeros(numExperiments, totalNumQueries);
% For storing cumulative regret values
boKDCumRegrets = zeros(numExperiments, totalNumQueries);
boAddCumRegrets = zeros(numExperiments, totalNumQueries, numdCands);
boEICumRegrets = zeros(numExperiments, totalNumQueries, numdCands);
randCumRegrets = zeros(numExperiments, totalNumQueries);

% First Call Direct
fprintf('First Running DiRect.\m');
diRectOpts.maxevals = totalNumQueries;
diRectOpts.maxits = inf;
diRectOpts.showits = true;
[~, ~, diRectHist, ~, diRectHistory] = diRectWrap(func, bounds, diRectOpts);
[diRectSimpleRegret, diRectCumRegret] = getRegrets(trueMaxVal, diRectHistory);
fprintf('DiRect OptVal: %0.4f\n', diRectHist(end));

for expIter = 1:numExperiments

  fprintf('\n==============================================================\n');
  fprintf('Experiment %d/ %d\nMaxVal: %0.4f\n', ...
    expIter, numExperiments, trueMaxVal);
  fprintf('Num DiRectEvals: %d\n', numDiRectEvals);
  fprintf('==============================================================\n');

  % Known true decomposition
  fprintf('\nKnown Decomposition\n');
  boKDParams.noises = 0 * ones(trueNumGroups, 1);
  [~, ~, ~, valHistKD] = ...
    addGPBO(func, trueDecomp, bounds, numIters, boKDParams);
  [sR, cR] = getRegrets(trueMaxVal, valHistKD);
  boKDHistories(expIter, :) = valHistKD';
  boKDSimpleRegrets(expIter, :) = sR';
  boKDCumRegrets(expIter, :) = cR';

  % For the candidates in numDimsPerGroupCands
  for candIter = 1:numdCands
    fprintf('\nUsing an arbitrary %d/ %d decomposition\n', ...
      numDimsPerGroupCands(candIter), numDims );
    [decompAdd, boAddParamsCurr, numCurrGroups] = ...
      preprocessDecomposition(numDims, numDimsPerGroupCands(candIter), ...
                              boAddParams, true);
    boAddParamsCurr.diRectParams.maxevals = ...
                                ceil(0.9 * numDiRectEvals/numCurrGroups);
    [~, ~, ~, valHistAdd] = ...
      addGPBO(func, decompAdd, bounds, numIters, boAddParamsCurr);
    [sR, cR] = getRegrets(trueMaxVal, valHistAdd);
    boAddHistories(expIter, :, candIter) = valHistAdd';
    boAddSimpleRegrets(expIter, :, candIter) = sR';
    boAddCumRegrets(expIter, :, candIter) = cR';
  end

  % GP - EI
  fprintf('\nGP-EI\n');
  [~, ~, ~, valHistEI] = ...
    addGPBO(func, boEIDecomp, bounds, numIters, boEIParams);
  [sR, cR] = getRegrets(trueMaxVal, valHistEI);
  boEIHistories(expIter, :, candIter) = valHistEI;
  boEISimpleRegrets(expIter, :, candIter) = sR';
  boEICumRegrets(expIter, :, candIter) = cR';

  % Random
  randQueries = bsxfun(@plus, ...
    bsxfun(@times, rand(totalNumQueries, numDims), ...
      (bounds(:,2) - bounds(:,1))' ), bounds(:,1)' );
  randHistories(expIter, :) = func(randQueries)';
  [sR, cR] = getRegrets(trueMaxVal, randHistories(expIter, :)');
  randSimpleRegrets(expIter, :) = sR';
  randCumRegrets(expIter, :) = cR';
  
end

% Plot the results
plotGPBResults;

