function [nextPt, nextPtMean, nextPtStd, nextPtUtil] = ...
  boGetNextQueryPt(params, gpFuncH, boVals, bounds)
% This is what this function should do. It should maximize the utility 

  % First create the utility function to be maximized
  if strcmp(params.utilityFunc, 'EI')
    utility = @(t) getEIUtility(t, gpFuncH, max(boVals));
  else
    utility = @(t) params.utilityFunc(t, gpFuncH);
  end

  % Now run DiRect
  [nextPtUtil, nextPt, hist] = diRectWrap(utility, bounds, params.diRectParams);

  % Now get the mean and variance of nextPt
  [nextPtMean, nextPtStd] = gpFuncH(nextPt);

end
