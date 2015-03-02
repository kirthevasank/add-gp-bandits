function K = augKernelNoise(X, coords, bw, scale, noise)
% Same as above, but adds the noise too
  K = augKernel(X, X, coords, bw, scale) + noise * eye(size(X,1));
end
