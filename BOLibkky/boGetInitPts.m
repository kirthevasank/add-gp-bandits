function initPts = boGetInitPts(bounds, numInitPts)
% This is just a utility function to obtain a set of initial points. We will be
% using it across all our BO implementations.

  useKMPP = false;
  numDims = size(bounds, 1);

  if numInitPts == 1
    initPts = [ ( bounds(:,2) + bounds(:,1) )/2 ]';
  else

    if useKMPP
      if params.numInitPts < 10, numRandPts = 200;
      else numRandPts = 10*numInitPts;
      end

      initPtCandidates = bsxfun(@plus, ...
        bsxfun(@times, rand(numRandPts, numDims), ...
               (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
      [~, initPts] = kmeanspp(initPtCandidates, params.numInitPts);
    else
      initPts = bsxfun(@plus, ...
        bsxfun(@times, rand(numInitPts, numDims), ...
               (bounds(:,2) - bounds(:,1))' ), bounds(:,1)');
    end
  end

end
