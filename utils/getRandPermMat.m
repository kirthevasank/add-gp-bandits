function P = getRandPermMat(D)
  P = zeros(D);
  shuffleOrder = randperm(D);
  for i = 1:D
    P(i, shuffleOrder(i)) = 1;
  end
end
