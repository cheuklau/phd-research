function CheckMoments(quadFrom, solFrom, quadTo, solTo)

% 0th and 1st moments from quadrature you are mapping from
numDirsFrom = size(quadFrom, 1);
zerothFrom = 0;
firstXfrom = 0;
firstYfrom = 0;
firstZfrom = 0;
secondx2y2From = 0;
secondz2From = 0;
secondxyFrom = 0;
secondxzFrom = 0;
secondyzFrom = 0;
for i = 1 : numDirsFrom
      x = quadFrom{i}(1);
      y = quadFrom{i}(2);
      z = quadFrom{i}(3);
      wt = quadFrom{i}(4);      
      sol = solFrom(i);
      zerothFrom = zerothFrom + wt * sol;
      firstXfrom = firstXfrom + wt * x * sol;
      firstYfrom = firstYfrom + wt * y * sol;
      firstZfrom = firstZfrom + wt * z * sol;
      secondx2y2From = secondx2y2From + wt * (x ^ 2 - y ^ 2) * sol;
      secondz2From = secondz2From + wt * z ^ 2 * sol;
      secondxyFrom = secondxyFrom + wt * x * y * sol;
      secondxzFrom = secondxzFrom + wt * x * z * sol;
      secondyzFrom = secondyzFrom + wt * y * z * sol;
end

% Zeroth and first moments from quadrature you are mapping to
numDirsTo = size(quadTo, 1);
numSubSqTo = numDirsTo / (9 * 3);
zerothTo = 0;
firstXto = 0;
firstYto = 0;
firstZto = 0;
secondx2y2To = 0;
secondz2To = 0;
secondxyTo = 0;
secondxzTo = 0;
secondyzTo = 0;
counter = 1;
for i = 1 : 3
    for j = 1 : numSubSqTo
        for k = 1 : 9
            x = quadTo{counter}(1);
            y = quadTo{counter}(2);
            z = quadTo{counter}(3);
            wt = quadTo{counter}(4);
            sol = solTo{i}{j}(k);
            zerothTo = zerothTo + wt * sol;
            firstXto = firstXto + wt * x * sol;
            firstYto = firstYto + wt * y * sol;
            firstZto = firstZto + wt * z * sol;
            secondx2y2To = secondx2y2To + wt * (x ^ 2 - y ^ 2) * sol;
            secondz2To = secondz2To + wt * z ^ 2 * sol;
            secondxyTo = secondxyTo + wt * x * y * sol;
            secondxzTo = secondxzTo + wt * x * z * sol;
            secondyzTo = secondyzTo + wt * y * z * sol;
            counter = counter + 1;
        end
    end
end

% Calculate relative different between the two moments
errZeroth = abs((zerothFrom - zerothTo) / zerothFrom);
errFirstX = abs((firstXfrom - firstXto) / firstXfrom);
errFirstY = abs((firstYfrom - firstYto) / firstYfrom);
errFirstZ = abs((firstZfrom - firstZto) / firstZfrom);
if abs(secondx2y2From) > 1e-6 && abs(secondx2y2To) > 1e-6
    errSecondx2y2 = abs((secondx2y2From - secondx2y2To) / secondx2y2From);
else
    errSecondx2y2 = 0;
end
errSecondz2 = abs((secondz2From - secondz2To) / secondz2From);
errSecondxy = abs((secondxyFrom - secondxyTo) / secondxyFrom);
errSecondxz = abs((secondxzFrom - secondxzTo) / secondxzFrom);
errSecondyz = abs((secondyzFrom - secondyzTo) / secondyzFrom);

% Print error
fprintf('Error in zeroth moment     : %e \n', errZeroth);
fprintf('Error in first moment (mu) : %e \n', errFirstX);
fprintf('Error in first moment (eta): %e \n', errFirstY);
fprintf('Error in first moment (xi) : %e \n', errFirstZ);
fprintf('Error in second moment (mu ^ 2 - eta ^ 2): % e \n', errSecondx2y2);
fprintf('Error in second moment (xi ^ 2): % e \n', errSecondz2);
fprintf('Error in second moment (mu * eta): % e \n', errSecondxy);
fprintf('Error in second moment (mu * xi): % e \n', errSecondxz);
fprintf('Error in second moment (eta * xi): % e \n', errSecondyz);
fprintf('RMS error in first and second moments: %f \n', 100*sqrt((errFirstX ^ 2 + ...
    errFirstY ^ 2 + errFirstZ ^ 2 + errSecondx2y2 ^ 2 + errSecondxy ^ 2 + ...
    errSecondxz ^ 2 + errSecondyz ^ 2 + errSecondz2 ^ 2) / 8));

end