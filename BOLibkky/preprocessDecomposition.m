function [decomp, params, numGroups] = ...
  preprocessDecomposition(numDims, numDimsPerGroup, params, addRemainingDims)

  % First set the Decomposition
  if ~isfield(params, 'decompStrategy') | isempty(params.decompStrategy)
    params.decompStrategy = 'partialLearn';
  end

  if ~exist('addRemainingDims', 'var')
    addRemainingDims = true;
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
    % This is full (naive) BO
    params.decompStrategy = 'known';
    decomp = cell(1,1);
    decomp{1} = 1:numDims;
    params.noises = 0 * ones(numGroups, 1);

  elseif strcmp(params.decompStrategy, 'known')
    decomp = cell(numGroups, 1);
    params.noises = 0 * ones(numGroups, 1);
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

