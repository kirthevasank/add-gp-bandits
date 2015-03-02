function [decomp, boAddParams, numGroups] = ...
  getDecompForParams(numDims, numDimsPerGroup, boAddParams, addRemainingDims)

  if ~exist('addRemainingDims', 'var')
    addRemainingDims = false;
  end

  if addRemainingDims
    numGroups = ceil(numDims/numDimsPerGroup);
    numRemDims = numDims - numDimsPerGroup * (numGroups-1);
    if numRemDims == 0
      addRemainingDims = false;
    end
  else
    numGroups = floor(numDims/numDimsPerGroup);
  end

  if numDimsPerGroup == numDims
    % This is the full BO
    boAddParams.decompStrategy = 'known';
    decomp = cell(1,1);
    decomp{1} = 1:numDims;
    boAddParams.noises = 0 * ones(numGroups, 1);

  elseif strcmp(boAddParams.decompStrategy, 'known')
    decomp = cell(numGroups, 1);
    boAddParams.noises = 0 * ones(numGroups, 1);
    for i = 1:numGroups
      decomp{i} = ( (i-1)*numDimsPerGroup+1 : min(i*numDimsPerGroup, numDims) );
    end

  elseif addRemainingDims
    decomp = [ numDimsPerGroup * ones(numGroups-1, 1); numRemDims];
  
  else
    decomp.d = numDimsPerGroup;
    decomp.M = numGroups;

  end

end

