function [hist, simpleReg, cumReg] = ...
  getDiRectResults(diRectHist, trueMaxVal, totalNumQueries)
% This is just a wrapper function to get the results for DiRect

  maxValKnots = [0; round(diRectHist(:,2))];
  numKnots = size(diRectHist, 1);
  numQueries = maxValKnots(end);

  hist = zeros(numQueries, 1);
  for i = 1: numKnots
    idxs = (maxValKnots(i)+1):maxValKnots(i+1);
    hist(idxs) = diRectHist(i, 3);
  end

  % Finally truncate at totalNumQueries
  hist = hist(1:totalNumQueries);
  [simpleReg, cumReg] = getRegrets(trueMaxVal, hist);

end

