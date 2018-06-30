% Generate QDFE-ratio quadrature
function quadrature = QDFEsa(squareInfo)

% Retrieve needed square properties
minSubSubX    = squareInfo.minSubSubX;
maxSubSubX    = squareInfo.maxSubSubX;
minSubSubY    = squareInfo.minSubSubY;
maxSubSubY    = squareInfo.maxSubSubY;
subSubLengthX = squareInfo.subSubLengthX;
subSubLengthY = squareInfo.subSubLengthY;
surfaceArea   = squareInfo.surfaceArea;
integrations  = squareInfo.integrations;
numSubSq      = squareInfo.numSubSq;

% Initialize storage
xPos         = cell(numSubSq, 1);
yPos         = cell(numSubSq, 1);
gamma        = cell(numSubSq, 3);
theta        = cell(numSubSq, 3);
weights      = cell(numSubSq, 1);
basis        = cell(numSubSq, 1);
basis2       = cell(numSubSq, 1);
basis3       = cell(numSubSq, 1);
weightMat    = zeros(8, 8);
tol_store    = zeros(numSubSq, 1);
factor_store = zeros(numSubSq, 1);

% Initial ratio for QDFE-ratio method
% rhoFixed is the non-corner sub-square ratio
[rho_corner_all, rho_side_all] = InitialRatio(squareInfo);

% Begin secant method
maxIter = 100;

% Go through each sub-square
for iSub = 1 : numSubSq
        
    % Surface area without center sub-sub-square
    surfaceAreaTemp(1) = surfaceArea{iSub}(1);
    surfaceAreaTemp(2) = surfaceArea{iSub}(2);
    surfaceAreaTemp(3) = surfaceArea{iSub}(3);
    surfaceAreaTemp(4) = surfaceArea{iSub}(4);
    surfaceAreaTemp(5) = surfaceArea{iSub}(6);
    surfaceAreaTemp(6) = surfaceArea{iSub}(7);
    surfaceAreaTemp(7) = surfaceArea{iSub}(8);
    surfaceAreaTemp(8) = surfaceArea{iSub}(9);
    
    % Initial ratio for current sub-square
    rho_corner = rho_corner_all(iSub);
    rho_side   = rho_side_all(iSub);
    
    % Reset perturbation factor, delta tolerance and outer convergence
    factor = 1e-1;
    tol    = 1e-8;
    convergedOuter = 0;
    
    % Outer iteration
    while convergedOuter == 0 && factor >= 1e-12 && tol <= 1e-4
        
        % Initial ratio
        rhoPrev = [...
            rho_corner + factor, rho_side - factor, rho_corner + factor, ...
            rho_side   - factor, rho_side - factor, ...
            rho_corner + factor, rho_side - factor, rho_corner + factor];
        
        % Next ratio guess
        rhoNext = [...
            rho_corner - factor, rho_side + factor, rho_corner - factor, ...
            rho_side   + factor, rho_side + factor, ...
            rho_corner - factor, rho_side + factor, rho_corner - factor];
        
        % Sub-sub-square length
        
        % Reset inner iteartion counter and convergence
        counter        = 0;
        convergedInner = 0;
        
        % Inner iteration for current sub-square
        while counter < maxIter && convergedInner == 0
            
            % Initialize converged tracker
            convergedID = [0, 0, 0, 0, 0, 0, 0, 0];
            
            % Quadrature location and weights using initial ratio
            xPosStore = [...
                maxSubSubX{iSub}(1) - subSubLengthX{iSub}(1) * rhoPrev(1), ...
                (maxSubSubX{iSub}(2) + minSubSubX{iSub}(2)) / 2, ...
                minSubSubX{iSub}(3) + subSubLengthX{iSub}(3) * rhoPrev(3), ...
                ...
                maxSubSubX{iSub}(4) - subSubLengthX{iSub}(4) * rhoPrev(4), ...
                (maxSubSubX{iSub}(5) + minSubSubX{iSub}(5)) / 2, ...
                minSubSubX{iSub}(6) + subSubLengthX{iSub}(6) * rhoPrev(5), ...
                ...
                maxSubSubX{iSub}(7) - subSubLengthX{iSub}(7) * rhoPrev(6), ...
                (maxSubSubX{iSub}(8) + minSubSubX{iSub}(8)) / 2, ...
                minSubSubX{iSub}(9) + subSubLengthX{iSub}(9) * rhoPrev(8)];
            
            yPosStore = [...
                maxSubSubY{iSub}(1) - subSubLengthY{iSub}(1) * rhoPrev(1), ...
                maxSubSubY{iSub}(2) - subSubLengthY{iSub}(2) * rhoPrev(2), ...
                maxSubSubY{iSub}(3) - subSubLengthY{iSub}(3) * rhoPrev(3), ...
                ...
                (maxSubSubY{iSub}(4) + minSubSubY{iSub}(4)) / 2, ...
                (maxSubSubY{iSub}(5) + minSubSubY{iSub}(5)) / 2, ...
                (maxSubSubY{iSub}(6) + minSubSubY{iSub}(6)) / 2, ...
                ...
                minSubSubY{iSub}(7) + subSubLengthY{iSub}(7) * rhoPrev(6), ...
                minSubSubY{iSub}(8) + subSubLengthY{iSub}(8) * rhoPrev(7), ...
                minSubSubY{iSub}(9) + subSubLengthY{iSub}(9) * rhoPrev(8)];
            [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
            [constants, constants2, constants3] = Basis(gammaStore, thetaStore);
            if constants == zeros(9)
                convergedInner = 1;
            end
            weightsStore = Weights(integrations{iSub}, constants);
            
            % Store weights in weight matrix
            weightMat(1, 1) = weightsStore(1);
            weightMat(1, 2) = weightsStore(2);
            weightMat(1, 3) = weightsStore(3);
            weightMat(1, 4) = weightsStore(4);
            weightMat(1, 5) = weightsStore(6);
            weightMat(1, 6) = weightsStore(7);
            weightMat(1, 7) = weightsStore(8);
            weightMat(1, 8) = weightsStore(9);
            
            % Form secant matrix moving one rho at a time
            numConv = 0;
            for m = 1 : 8
                rhoTest = rhoPrev;
                rhoTest(m) = rhoNext(m);
                xPosTemp = [...
                    maxSubSubX{iSub}(1) - subSubLengthX{iSub}(1) * rhoTest(1), ...
                    (maxSubSubX{iSub}(2) + minSubSubX{iSub}(2)) / 2, ...
                    minSubSubX{iSub}(3) + subSubLengthX{iSub}(3) * rhoTest(3), ...
                    ...
                    maxSubSubX{iSub}(4) - subSubLengthX{iSub}(4) * rhoTest(4), ...
                    (maxSubSubX{iSub}(5) + minSubSubX{iSub}(5)) / 2, ...
                    minSubSubX{iSub}(6) + subSubLengthX{iSub}(6) * rhoTest(5), ...
                    ...
                    maxSubSubX{iSub}(7) - subSubLengthX{iSub}(7) * rhoTest(6), ...
                    (maxSubSubX{iSub}(8) + minSubSubX{iSub}(8)) / 2, ...
                    minSubSubX{iSub}(9) + subSubLengthX{iSub}(9) *  rhoTest(8)];
                
                yPosTemp = [...
                    maxSubSubY{iSub}(1) - subSubLengthY{iSub}(1) * rhoTest(1), ...
                    maxSubSubY{iSub}(2) - subSubLengthY{iSub}(2) *  rhoTest(2), ...
                    maxSubSubY{iSub}(3) - subSubLengthY{iSub}(3) *  rhoTest(3), ...
                    ...
                    (maxSubSubY{iSub}(4) + minSubSubY{iSub}(4)) / 2, ...
                    (maxSubSubY{iSub}(5) + minSubSubY{iSub}(5)) / 2, ...
                    (maxSubSubY{iSub}(6) + minSubSubY{iSub}(6)) / 2, ...
                    ...
                    minSubSubY{iSub}(7) + subSubLengthY{iSub}(7) *  rhoTest(6), ...
                    minSubSubY{iSub}(8) + subSubLengthY{iSub}(8) *  rhoTest(7), ...
                    minSubSubY{iSub}(9) + subSubLengthY{iSub}(9) *  rhoTest(8)];
                [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                [constants, ~, ~] = Basis(gammaTemp, thetaTemp);
                if constants == zeros(9)
                    convergedInner = 1;
                end
                weightsTemp = Weights(integrations{iSub}, constants);
                weightMat(m + 1, 1) = weightsTemp(1);
                weightMat(m + 1, 2) = weightsTemp(2);
                weightMat(m + 1, 3) = weightsTemp(3);
                weightMat(m + 1, 4) = weightsTemp(4);
                weightMat(m + 1, 5) = weightsTemp(6);
                weightMat(m + 1, 6) = weightsTemp(7);
                weightMat(m + 1, 7) = weightsTemp(8);
                weightMat(m + 1, 8) = weightsTemp(9);
                
                % Check if any ratios are super-converged
                if weightMat(m + 1, 1) == weightMat(1, 1) &&...
                        weightMat(m + 1, 2) == weightMat(1, 2) &&...
                        weightMat(m + 1, 3) == weightMat(1, 3) &&...
                        weightMat(m + 1, 4) == weightMat(1, 4) &&...
                        weightMat(m + 1, 5) == weightMat(1, 5) &&...
                        weightMat(m + 1, 6) == weightMat(1, 6) &&...
                        weightMat(m + 1, 7) == weightMat(1, 7) &&...
                        weightMat(m + 1, 8) == weightMat(1, 8)
                    numConv = numConv + 1;
                    convergedID(m) = 1;
                end
            end
            
            % No points are converged
            delta = [0, 0, 0, 0, 0, 0, 0, 0];
            if numConv == 0
                deriv = zeros(8);
                difference = zeros(1, 8);
                for i = 1 : 8
                    for j = 1 : 8
                        deriv(i, j) = (weightMat(j + 1, i) - weightMat(1, i)) ...
                            / (rhoNext(j) - rhoPrev(j));
                    end
                    difference(i) = surfaceAreaTemp(i) - weightMat(1, i);
                end
                if rcond(deriv) < eps
                    convergedInner = 1;
                else
                    delta = deriv \ difference';
                end
                
                % If one through seven points are converged
            elseif numConv > 0 && numConv < 7
                numUnconv = 8 - numConv;
                deriv = zeros(numUnconv);
                difference = zeros(1, numUnconv);
                indexI = 1;
                for i = 1 : 8
                    if convergedID(i) ~= 1
                        indexJ = 1;
                        for j = 1 : 8
                            if convergedID(j) ~= 1
                                deriv(indexI, indexJ) = ...
                                    (weightMat(j + 1, i) - weightMat(1, i)) ...
                                    / (rhoNext(j) - rhoPrev(j));
                                indexJ = indexJ + 1;
                            end
                        end
                        difference(indexI) = surfaceAreaTemp(i) - weightMat(1, i);
                        indexI = indexI + 1;
                    end
                end
                if rcond(deriv) < eps
                    convergedInner = 1;
                else
                    deltaTemp = deriv \ difference';
                    indexTemp = 1;
                    for i = 1 : 8
                        if convergedID(i) ~= 1
                            delta(i) = deltaTemp(indexTemp);
                            indexTemp = indexTemp + 1;
                        end
                    end
                end
                
                % Seven points are converged
            elseif numConv == 7
                for i = 1 : 8
                    if convergedID(i) == 0;
                        unconvergedID = i;
                    end
                end
                deriv = (weightMat(unconvergedID + 1, unconvergedID) - ...
                    weightMat(1, unconvergedID)) / (rhoNext(unconvergedID) - ...
                    rhoPrev(unconvergedID));
                difference = surfaceAreaTemp(unconvergedID) - weightMat(1, unconvergedID);
                delta(unconvergedID) = difference / deriv;
                
                % All point are converged
            elseif numConv == 8
                delta = [0, 0, 0, 0, 0, 0, 0, 0];
                
                % Exception
            else
                disp('Error in numConv definition! \n')
            end
            
            % Next ratio guess
            rhoTemp = zeros(1, 8);
            for i = 1 : 8
                rhoTemp(i) = rhoPrev(i) + delta(i);
            end
            
            % Check for convergence
            %{
            if abs(surfaceAreaTemp(1) - weightMat(1, 1)) / surfaceAreaTemp(1) < tol && ...
               abs(surfaceAreaTemp(2) - weightMat(1, 2)) / surfaceAreaTemp(2) < tol && ...
               abs(surfaceAreaTemp(3) - weightMat(1, 3)) / surfaceAreaTemp(3) < tol && ...
               abs(surfaceAreaTemp(4) - weightMat(1, 4)) / surfaceAreaTemp(4) < tol && ...
               abs(surfaceAreaTemp(5) - weightMat(1, 5)) / surfaceAreaTemp(5) < tol && ...
               abs(surfaceAreaTemp(6) - weightMat(1, 6)) / surfaceAreaTemp(6) < tol && ...
               abs(surfaceAreaTemp(7) - weightMat(1, 7)) / surfaceAreaTemp(7) < tol && ...
               abs(surfaceAreaTemp(8) - weightMat(1, 8)) / surfaceAreaTemp(8) < tol
            %}
            if abs(delta(1) / rhoPrev(1)) < tol && ...
                    abs(delta(2) / rhoPrev(2)) < tol && ...
                    abs(delta(3) / rhoPrev(3)) < tol && ...
                    abs(delta(4) / rhoPrev(4)) < tol && ...
                    abs(delta(5) / rhoPrev(5)) < tol && ...
                    abs(delta(6) / rhoPrev(6)) < tol && ...
                    abs(delta(7) / rhoPrev(7)) < tol && ...
                    abs(delta(8) / rhoPrev(8)) < tol && ...
                    convergedInner == 0
                convergedInner = 1;
                convergedOuter = 1;
            else
                rhoPrev = rhoNext;
                for i = 1 : 8
                    rhoNext(i) = rhoTemp(i);
                end
            end
            
            counter = counter + 1;
            
        end
        
        % If factor is minimized then increase tolerance
        if abs(factor - 1e-12) < eps && tol < 1e-5
            tol = tol * 10;
            factor = 1e-1;
        else
            % Update factor
            factor = factor / 10;
        end
        
    end
    
    % Store factor of each sub-square
    factor_store(iSub) = factor * 10;
    
    % Store error of each sub-square
    tol_store(iSub) = tol;
    
    % Check if current sub-square failed
    if convergedOuter == 0;
        error('QDFE-sa failed!');
    end    
    
    % Store quadrature
    xPos{iSub}     = xPosStore;
    yPos{iSub}     = yPosStore;
    gamma{iSub, 1} = gammaStore;
    theta{iSub, 1} = thetaStore;
    weights{iSub}  = weightsStore;
    basis{iSub}    = constants;
    basis2{iSub}   = constants2;
    basis3{iSub}   = constants3;
    
end

[min_val, min_pos] = min(factor_store);
fprintf('The minimum factor is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(factor_store);
fprintf('The maximum factor is: %E located in sub-square: %i \n', max_val, max_pos);

[min_val, min_pos] = min(tol_store);
fprintf('The minimum tol is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(tol_store);
fprintf('The maximum tol is: %E located in sub-square: %i \n', max_val, max_pos);

% Normalize weights to 4pi
tot_weight = 0;
counter = 1;
for i = 1 : numSubSq
    for j = 1 : 9
        tot_weight = tot_weight + weights{i}(j);
        counter = counter + 1;
    end
end
norm = (4 * pi / 24) / tot_weight;
for i = 1 : numSubSq
    for j = 1 : 9
        weights{i}(j) = weights{i}(j) * norm;
    end
end

% Maximum and average surface area error
test = zeros(numSubSq, 1);
counter = 1;
for i = 1 : numSubSq
    for j = 1 : 9
        test(counter) = abs(weights{i}(j) - surfaceArea{i}(j)) / surfaceArea{i}(j);
        counter = counter + 1;
    end
end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

% Rotate to other three faces
for i = 1 : numSubSq
    for j = 1 : 9
        gamma{i, 2}(j) = (pi / 2) - atan(sqrt(3) * yPos{i}(j));
        theta{i, 2}(j) = (pi / 2) - atan(sqrt(3) * xPos{i}(j) /...
            sqrt(1 + 3 * yPos{i}(j) ^ 2));
        gamma{i, 3}(j) = atan(yPos{i}(j) / xPos{i}(j));
        theta{i, 3}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos{i}(j) ^ 2 + ...
            xPos{i}(j) ^ 2)));
    end
end

% Store QDFE-sa quadrature data
quadrature.gamma   = gamma;
quadrature.theta   = theta;
quadrature.weights = weights;
quadrature.basis   = basis;
quadrature.basis2  = basis2;
quadrature.basis3  = basis3;

end