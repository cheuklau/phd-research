function CheckMoments(quadFrom, solFrom, quadTo, solTo)

% 0th and 1st moments from quadrature you are mapping from
weightFrom = quadFrom.weight;
muFrom     = quadFrom.mu;
zerothMomentFrom = 0;
firstMomentFrom  = 0;
numRegionsFrom = size(solFrom, 1);
for i = 1 : numRegionsFrom
   for j = 1 : 2
      zerothMomentFrom = zerothMomentFrom + weightFrom{i}(j) * ...
          solFrom{i}(j);
      firstMomentFrom = firstMomentFrom + weightFrom{i}(j) * ...
          muFrom{i}(j) * solFrom{i}(j);       
   end    
end

% Zeroth and first moments from quadrature you are mapping to
weightTo = quadTo.weight;
muTo     = quadTo.mu;
zerothMomentTo = 0;
firstMomentTo = 0;
numRegionsTo = size(solTo, 1);
for i = 1 : numRegionsTo
   for j = 1 : 2
      zerothMomentTo = zerothMomentTo + weightTo{i}(j) * ...
          solTo{i}(j);
      firstMomentTo = firstMomentTo + weightTo{i}(j) * ...
          muTo{i}(j) * solTo{i}(j);
   end    
end

% Calculate relative different between the two moments
errZeroth = abs((zerothMomentFrom - zerothMomentTo) / zerothMomentFrom);
errFirst = abs((firstMomentFrom - firstMomentTo) / firstMomentFrom);
fprintf('error in zeroth moment: %e \n', errZeroth);
fprintf('error in first moment: %e \n', errFirst);
end