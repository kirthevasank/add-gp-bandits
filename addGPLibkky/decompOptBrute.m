function [A, optVal] = decompOptBrute(func, D, d, M)

  % This searches over all permutations
  if D > 12
    error('D is too large');
  end

  p = d*M;
  numCombs = factorial(D)/ ( factorial(M) * factorial(d)^M * factorial(D-p) );
%   nc2 = exp(gammaln(D+1) - M*gammaln(d+1) - gammaln(D-p+1) - gammaln(M+1));
%   fprintf('Num Combinations: %d, %f,\n', numCombs, nc2);
%   pause,
  
  optVal = +inf;
  counts = zeros(M);
  combCounter = 0;

  relevantCombs = nchoosek(1:D, p);
  numRelCombs = size(relevantCombs, 1);
  count = 0;
  for relCoordIter = 1:numRelCombs

    var1 = relevantCombs(relCoordIter, :);

    % 1st group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g1 = nextChooseTerminate(p, d);
    while true
      aa1 = g1(); a1 = var1(aa1);
      if isSizeZero(a1), break, end
      if M >= 2, var2 = setdiff(var1, a1);
      else, var2 = []; end

      % 2nd group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      g2 = nextChooseTerminate(p-d, d);
      while true 
        aa2 = g2(); a2 = var2(aa2);
        if isSizeZero(a2), break, end
        if M>= 3, var3 = setdiff(var2, a2);
        else, var3 = []; end
        
        % 3rd group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        g3 = nextChooseTerminate(p-2*d, d);
        while true
          aa3 = g3(); a3 = var3(aa3);
          if isSizeZero(a3), break, end
          if M>=4, var4 = setdiff(var3, a3);
          else, var4 = []; end

          % 4th group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          g4 = nextChooseTerminate(p-3*d, d);
          while true
            aa4 = g4(); a4 = var4(aa4);
            if isSizeZero(a4), break, end
            
            b = [a1 a2 a3 a4];
            candA = getA(b, D, p);
            val = func(candA);
            if val < optVal
              A = candA;
              optVal = val;
            end
            count = count + 1;
%             fprintf('# %d, %s, %.4f\n', count, mat2str(b), val);
%             candA,
          end
          % 4th group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        % 3rd group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
      % 2nd group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    % 1st group %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end

end


function A = getA(order, D, p)
  A = zeros(D, p);
  for i = 1:p
    A(order(i), i) = 1;
  end
end


function t = isSizeZero(A)
  t = (size(A,1) == 0) & (size(A,2) == 0);
end

