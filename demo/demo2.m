% Experiment Set up for Bayesian Optimization and GP Bandits

close all;
clear all;
addpath ../addGPLibkky/
addpath ../BOLibkky/
addpath ../utils/
warning off;
rng('shuffle');

% Problem Set up
% --------------
numExperiments = 2; % Running 2 experiments to generate the error bars.
numDims = 36; trueNumDimsPerGroup = 11;
% numDims = 40; trueNumDimsPerGroup = 18;
% Get the function
[func, funcProperties] = getAdditiveFunction(numDims, trueNumDimsPerGroup);
bounds = funcProperties.bounds;
trueMaxVal = funcProperties.maxVal;

% EXPERIMENT PARAMETERS: At minimum, you need to set these two values.
% --------------------------------------------------------------------
d = 5; % Add-GP-UCB with maximum group size d;
numIters = 1000; % The number of iterations for BO/Bandits.
params = struct();
params.decompStrategy = 'stoch1';

% Call additive GP Bandits
% ------------------------
% First call the following function to set up the hyper parameters
fprintf('Running BO (Maxval: %0.4f.)\n', funcProperties.maxVal)
[decomp, boAddParams] = ...
  preprocessDecomposition(numDims, d, params, true);
% Call adGPBO
[maxVal, maxPt, boQueries, boVals, history] = ...
  addGPBO(func, decomp, bounds, numIters, boAddParams);

% Obtain the simple and cumulative regrets and plot them.
[sR, cR] = getRegrets(trueMaxVal, history);
figure; plot(sR); title('Simple Regrets');
figure; plot(cR); title('(Time Averaged) - Cumulative Regrets');

