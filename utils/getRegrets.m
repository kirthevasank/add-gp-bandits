function [simpleRegret, cumRegretByT, cumRegret] = getRegrets(maxVal, histories)
% If maxval is not finite, then we will just return the rewards

  numQueries = size(histories, 1);

  if isfinite(maxVal)
    diffs = maxVal - histories;
  else
    diffs = histories;
  end
  cumRegret = cumsum(diffs);
  cumRegretByT = cumRegret ./ (1:numQueries)';

  simpleRegret = zeros(numQueries, 1);
  simpleRegret(1) = diffs(1);
  for i = 2:numQueries
    if isfinite(maxVal)
      simpleRegret(i) = min( simpleRegret(i-1), diffs(i) );
    else
      simpleRegret(i) = max( simpleRegret(i-1), diffs(i) );
    end
  end

end

