function samples = sampleFromMultinomial(p, n)
  if ~exist('n', 'var')
    n = 1;
  end
  preSamples = mnrnd(1, p, n);
  [samples, ~] = find(preSamples');
end

