function As, probs = decompOptPartialStoch(func, D, dMax, M)

  numTrials = 5;
  logLikls = zeros(dMax, 1);
  As = cell(dMax, 1);

  for d = 1:dMax
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

    As{d} = A;
    dMax(d) = -currBestVal;

  end

  probs = exp(logLikls);
  probs = probs/sum(probs);

end

