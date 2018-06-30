function CheckMoments(quadFrom, solFrom, quadTo, solTo)

% 0th and 1st moments from quadrature you are mapping from
numDirsFrom = size(quadFrom, 1);
zerothFrom = 0;
firstXfrom = 0;
firstYfrom = 0;
firstZfrom = 0;
for i = 1 : numDirsFrom
      zerothFrom = zerothFrom + quadFrom{i}(4) * solFrom(i);
      firstXfrom = firstXfrom + quadFrom{i}(4) * quadFrom{i}(1) * solFrom(i);
      firstYfrom = firstYfrom + quadFrom{i}(4) * quadFrom{i}(2) * solFrom(i);
      firstZfrom = firstZfrom + quadFrom{i}(4) * quadFrom{i}(3) * solFrom(i);
end

% Zeroth and first moments from quadrature you are mapping to
numDirsTo = size(quadTo, 1);
numSubSqTo = numDirsTo / (4 * 3);
zerothTo = 0;
firstXto = 0;
firstYto = 0;
firstZto = 0;
counter = 1;
for i = 1 : 3
    for j = 1 : numSubSqTo
        for k = 1 : 4
            zerothTo = zerothTo + quadTo{counter}(4) * solTo{i}{j}(k);
            firstXto = firstXto + quadTo{counter}(4) * quadTo{counter}(1) * solTo{i}{j}(k);
            firstYto = firstYto + quadTo{counter}(4) * quadTo{counter}(2) * solTo{i}{j}(k);
            firstZto = firstZto + quadTo{counter}(4) * quadTo{counter}(3) * solTo{i}{j}(k);
            counter = counter + 1;
        end
    end
end

% Calculate relative different between the two moments
errZeroth = abs((zerothFrom - zerothTo) / zerothFrom);
errFirstX = abs((firstXfrom - firstXto) / firstXfrom);
errFirstY = abs((firstYfrom - firstYto) / firstYfrom);
errFirstZ = abs((firstZfrom - firstZto) / firstZfrom);
fprintf('Error in zeroth moment     : %e \n', errZeroth);
fprintf('Error in first moment (mu) : %e \n', errFirstX);
fprintf('Error in first moment (eta): %e \n', errFirstY);
fprintf('Error in first moment (xi) : %e \n', errFirstZ);


end