% Generate LDFE-sa quadrature
% Note: Convergence very sensitive to the first ratio guess.
% Some refinements fail depending on which sub-sub-square you chose
% to converge for the first ratio guess! Instead of a first ratio guess 
% based on LDFE-ratio we may want to just sweep the entire range.

function quadrature = LDFEsa(squareInfo, rhoInitAll)

% Retrieve needed square properties
midSubX      = squareInfo.midSubX;
midSubY      = squareInfo.midSubY;
minSubSubX   = squareInfo.minSubSubX;
maxSubSubX   = squareInfo.maxSubSubX;
minSubSubY   = squareInfo.minSubSubY;
maxSubSubY   = squareInfo.maxSubSubY;
surfaceArea  = squareInfo.surfaceArea;
integrations = squareInfo.integrations;
numSubSq     = squareInfo.numSubSq;

% Initialize storage
xPos      = cell(numSubSq, 1);
yPos      = cell(numSubSq, 1);
gamma     = cell(numSubSq, 3);
theta     = cell(numSubSq, 3);
weights   = cell(numSubSq, 1);
basis     = cell(numSubSq, 1);
weightMat = zeros(4, 4);
tol_store = zeros(numSubSq, 1);
factor_store = zeros(numSubSq, 1);

% Begin secant method
maxIter = 100;

% Go through each sub-square
for iSub = 1 : numSubSq
    
    fprintf('sub-square %i of %i \n', iSub, numSubSq);
    
    % Initial ratio for current sub-square
    rhoInit = rhoInitAll(iSub);
    
    % Initial factor resulting in a ratio of 1
    factorRestart = 1.0 - rhoInit;
    
    % Amount to change initial factor each iteration
    deltaFactorRestart = 1 / 1000;
    
    % Reinitialize outer iteration convergence tracker
    convergedOuter = 0;
    
    % Iterate until factor results in a ratio of -1
    while convergedOuter == 0 && factorRestart > -1 * rhoInit
        
        % Reset perturbation factor and tolerance
        factor = factorRestart;
        tol    = 1e-12;
        
        % Outer iteration
        while convergedOuter == 0 && factor > 1e-12 % tol < 1e-11
            
            % Initial ratio
            rhoPrev = [rhoInit, rhoInit, rhoInit];
            
            % Next ratio guess
            rhoNext = [rhoInit + factor, rhoInit + factor, rhoInit + factor];
            
            % Reset inner iteration counter and convergence
            counter = 0;
            convergedInner = 0;
            
            % Inner iteration for current sub-square
            while counter < maxIter && convergedInner == 0
                
                % Initialize converged tracker
                convergedID = [0, 0, 0];
                
                % Quadrature location and weights using previous ratio
                xPosStore = [...
                    midSubX(iSub) + (maxSubSubX{iSub}(1) - minSubSubX{iSub}(1)) * rhoPrev(1),...
                    midSubX(iSub) + (maxSubSubX{iSub}(2) - minSubSubX{iSub}(2)) * rhoPrev(2),...
                    midSubX(iSub) - (maxSubSubX{iSub}(3) - minSubSubX{iSub}(3)) * rhoPrev(3),...
                    midSubX(iSub) - (maxSubSubX{iSub}(4) - minSubSubX{iSub}(4)) * rhoInit];
                yPosStore = [...
                    midSubY(iSub) + (maxSubSubY{iSub}(1) - minSubSubY{iSub}(1)) * rhoPrev(1),...
                    midSubY(iSub) - (maxSubSubY{iSub}(2) - minSubSubY{iSub}(2)) * rhoPrev(2),...
                    midSubY(iSub) - (maxSubSubY{iSub}(3) - minSubSubY{iSub}(3)) * rhoPrev(3),...
                    midSubY(iSub) + (maxSubSubY{iSub}(4) - minSubSubY{iSub}(4)) * rhoInit];
                [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
                constants = Basis(gammaStore, thetaStore);
                if constants == zeros(4)
                    convergedInner = 1;
                end
                weightsStore = Weights(integrations{iSub}, constants);
                
                % Store weights in weight matrix
                for i = 1 : 4
                    weightMat(1, i) = weightsStore(i);
                end
                
                % Form secant matrix moving one ratio at a time
                numConv = 0;
                for m = 1 : 3
                    rhoTest = rhoPrev;
                    rhoTest(m) = rhoNext(m);
                    xPosTemp = [...
                        midSubX(iSub) + (maxSubSubX{iSub}(1) - minSubSubX{iSub}(1)) * rhoTest(1),...
                        midSubX(iSub) + (maxSubSubX{iSub}(2) - minSubSubX{iSub}(2)) * rhoTest(2),...
                        midSubX(iSub) - (maxSubSubX{iSub}(3) - minSubSubX{iSub}(3)) * rhoTest(3),...
                        midSubX(iSub) - (maxSubSubX{iSub}(4) - minSubSubX{iSub}(4)) * rhoInit];
                    yPosTemp = [...
                        midSubY(iSub) + (maxSubSubY{iSub}(1) - minSubSubY{iSub}(1)) * rhoTest(1),...
                        midSubY(iSub) - (maxSubSubY{iSub}(2) - minSubSubY{iSub}(2)) * rhoTest(2),...
                        midSubY(iSub) - (maxSubSubY{iSub}(3) - minSubSubY{iSub}(3)) * rhoTest(3),...
                        midSubY(iSub) + (maxSubSubY{iSub}(4) - minSubSubY{iSub}(4)) * rhoInit];
                    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                    constants = Basis(gammaTemp, thetaTemp);
                    if constants == zeros(4)
                        convergedInner = 1;
                    end
                    weightsTemp = Weights(integrations{iSub}, constants);
                    for j = 1 : 4
                        weightMat(m + 1, j) = weightsTemp(j);
                    end
                    
                    % Check if any ratios are super-converged
                    if weightMat(m + 1, 1) == weightMat(1, 1) &&...
                            weightMat(m + 1, 2) == weightMat(1, 2) &&...
                            weightMat(m + 1, 3) == weightMat(1, 3) &&...
                            weightMat(m + 1, 4) == weightMat(1, 4)
                        numConv = numConv + 1;
                        convergedID(m) = 1;
                    end
                end
                
                % No points are converged
                delta = [0, 0, 0];
                if numConv == 0
                    deriv = zeros(3);
                    difference = zeros(1,3);
                    for i = 1 : 3
                        for j = 1 : 3
                            deriv(i, j) = (weightMat(j + 1, i) - weightMat(1, i))...
                                / (rhoNext(j) - rhoPrev(j));
                        end
                        difference(i) = surfaceArea{iSub}(i) - weightMat(1, i);
                    end
                    if rcond(deriv) < eps
                        convergedInner = 1;
                    else
                        delta = deriv \ difference';
                    end
                    
                    % One point is converged
                elseif numConv == 1
                    deriv = zeros(2, 2);
                    difference = zeros(1, 2);
                    indexI = 1;
                    for i = 1 : 3
                        if convergedID(i) ~= 1
                            indexJ = 1;
                            for j = 1 : 3
                                if convergedID(j) ~= 1
                                    deriv(indexI, indexJ) =...
                                        (weightMat(j + 1, i) - weightMat(1, i))...
                                        / (rhoNext(j) - rhoPrev(j));
                                    indexJ = indexJ + 1;
                                end
                            end
                            difference(indexI) = surfaceArea{iSub}(i) - weightMat(1, i);
                            indexI = indexI + 1;
                        end
                    end
                    if rcond(deriv) < eps
                        convergedInner = 1;
                    else
                        deltaTemp = deriv \ difference';
                        indexTemp = 1;
                        for i = 1 : 3
                            if convergedID(i) ~= 1
                                delta(i) = deltaTemp(indexTemp);
                                indexTemp = indexTemp + 1;
                            end
                        end
                    end
                    
                    % Two points are converged
                elseif numConv == 2
                    for i = 1 : 3
                        if convergedID(i) == 0;
                            unconvergedID = i;
                        end
                    end
                    deriv = (weightMat(unconvergedID + 1, unconvergedID) - ...
                        weightMat(1, unconvergedID)) / (rhoNext(unconvergedID) - ...
                        rhoPrev(unconvergedID));
                    difference = surfaceArea{iSub}(unconvergedID) - weightMat(1, unconvergedID);
                    delta(unconvergedID) = difference / deriv;
                    
                    % All point are converged
                elseif numConv == 3
                    delta = [0, 0, 0];
                    
                    % Exception
                else
                    error('Error in numConv definition! \n');
                end
                
                % Next ratio guess
                rhoTemp = zeros(1, 3);
                for i = 1 : 3
                    rhoTemp(i) = rhoPrev(i) + delta(i);
                end
                
                % Check for inner convergence
                if abs(delta(1) / rhoPrev(1)) < tol && ...
                   abs(delta(2) / rhoPrev(2)) < tol && ...
                   abs(delta(3) / rhoPrev(3)) < tol && ...
                   convergedInner == 0
                    
                    convergedInner = 1;
                    
                    % Check for false convergence of ratios
                    if abs(surfaceArea{iSub}(1) - weightMat(1, 1)) / surfaceArea{iSub}(1) < 1e-5 && ...
                       abs(surfaceArea{iSub}(2) - weightMat(1, 2)) / surfaceArea{iSub}(2) < 1e-5 && ...
                       abs(surfaceArea{iSub}(3) - weightMat(1, 3)) / surfaceArea{iSub}(3) < 1e-5 && ...
                       abs(surfaceArea{iSub}(4) - weightMat(1, 4)) / surfaceArea{iSub}(4) < 1e-5 && ...
                       min(rhoPrev) > 0.0 && ...
                       max(rhoPrev) < 1.0
                        
                        convergedOuter = 1;
                        
                    end
                    
                    % Update next iteration ratios
                else
                    rhoPrev = rhoNext;
                    for i = 1 : 3
                        rhoNext(i) = rhoTemp(i);
                    end
                end
                
                % Update inner iteration counter
                counter = counter + 1;
                
            end
            
            % If factor is at minimum threshold
            if convergedOuter == 0
                                        
                factor = factor / 10;                                    
                
            end
            
        end
        
        % Decrease factor
        factorRestart = factorRestart - deltaFactorRestart;
        
        disp(factorRestart);
        
    end
    
    % Store factor of each sub-square
    factor_store(iSub) = factor;
    
    % Store error of each sub-square
    tol_store(iSub) = tol;
    
    % Check if current sub-square failed
    if convergedOuter == 0;
        error('LDFE-sa failed!');
    end
    
    % Store quadrature
    xPos{iSub}     = xPosStore;
    yPos{iSub}     = yPosStore;
    gamma{iSub, 1} = gammaStore;
    theta{iSub, 1} = thetaStore;
    weights{iSub}  = weightsStore;
    basis{iSub}    = constants;
    
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
    for j = 1 : 4
        tot_weight = tot_weight + weights{i}(j);
        counter = counter + 1;
    end
end
norm = (4 * pi / 24) / tot_weight;
for i = 1 : numSubSq
    for j = 1 : 4
        weights{i}(j) = weights{i}(j) * norm;
    end
end
fprintf('The 4pi normalization factor is: %E \n', norm);

% Actual surface area error
test = zeros(numSubSq * 4, 1);
tot_weight = 0;
counter = 1;
for i = 1 : numSubSq
    for j = 1 : 4
        tot_weight = tot_weight + weights{i}(j);
        test(counter) = abs(weights{i}(j) - surfaceArea{i}(j)) / surfaceArea{i}(j);
        counter = counter + 1;
    end
end
fprintf('the maximum surface area error is: %E \n', max(test));
fprintf('the average surface area error is: %E \n', mean(test));

% Rotate to other three faces
for i = 1 : numSubSq
    for j = 1 : 4
        gamma{i, 2}(j) = (pi / 2) - atan(sqrt(3) * yPos{i}(j));
        theta{i, 2}(j) = (pi / 2) - atan(sqrt(3) * xPos{i}(j) /...
            sqrt(1 + 3 * yPos{i}(j) ^ 2));
        gamma{i, 3}(j) = atan(yPos{i}(j) / xPos{i}(j));
        theta{i, 3}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos{i}(j) ^ 2 + ...
            xPos{i}(j) ^ 2)));
    end
end

% Store LDFE-sa quadrature data
quadrature.gamma   = gamma;
quadrature.theta   = theta;
quadrature.weights = weights;
quadrature.basis   = basis;

end