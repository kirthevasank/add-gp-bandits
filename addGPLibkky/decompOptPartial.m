function A = decompOptPartial(func, D, d, M)

  numTrials = 5;
  currBestVal = inf;

  for i = 1:numTrials
    P = getRandPermMat(D);
    val = func(P);
    if val < currBestVal
      A = P;
      currBestVal = val;
    end
  end

  if ~isfinite(currBestVal)
    A = getRandPermMat(D);
  end

end

