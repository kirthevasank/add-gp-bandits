function [C] = getHungarianCostMatrixOrth(A, costIdx, power)
% The objective is to take A - an orthogonal matrix and return the "closest"
% permutation matrix. We can treat this as an assignment problem where the cost
% of assigning the Ai to Zk (where Z is the output Perm matrix) is as follows:
% costIdx =1, Cik = 1 - |Akj|,  
% costIdx =2, Cik = 1 - |Akj| + \sum_{j=/=i} |Aij|
% costIdx =3, Cik = 1 - |Akj| + \sum_{j=/=i} |Aij| + \sum_{k=/=l} |Akl}

  % Prelims
  D = size(A, 1);
  p = size(A, 2);
  if ~exist('power', 'var')
    power = 2;
  end

  Aabs = abs(A).^power;
  C = 1 - Aabs;

  if costIdx >= 2
    for i = 1:D
      C1 = sum( Aabs(1:D ~= i, :) );
      C(i, :) = C(i,:) + C1;
    end

    if costIdx == 3
      for j = 1:p
        C2 = sum(Aabs(:, 1:p ~= j), 2);
        C(:, j) = C(:, j) + C2;
      end
    end
  end

end

