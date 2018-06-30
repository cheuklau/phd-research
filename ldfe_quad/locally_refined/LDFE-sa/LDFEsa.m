% Generate LDFE-center quadrature
function quadrature = LDFEsa(squareInfo, rhoInitAll1, rhoInitAll2, rhoInitAll3)

%% Retrieve needed square properties
midSubX1      = squareInfo.midSubX1;
midSubX2      = squareInfo.midSubX2;
midSubX3      = squareInfo.midSubX3;
midSubY1      = squareInfo.midSubY1;
midSubY2      = squareInfo.midSubY2;
midSubY3      = squareInfo.midSubY3;
minSubSubX1   = squareInfo.minSubSubX1;
minSubSubX2   = squareInfo.minSubSubX2;
minSubSubX3   = squareInfo.minSubSubX3;
maxSubSubX1   = squareInfo.maxSubSubX1;
maxSubSubX2   = squareInfo.maxSubSubX2;
maxSubSubX3   = squareInfo.maxSubSubX3;
minSubSubY1   = squareInfo.minSubSubY1;
minSubSubY2   = squareInfo.minSubSubY2;
minSubSubY3   = squareInfo.minSubSubY3;
maxSubSubY1   = squareInfo.maxSubSubY1;
maxSubSubY2   = squareInfo.maxSubSubY2;
maxSubSubY3   = squareInfo.maxSubSubY3;
integrations1 = squareInfo.integrations1;
integrations2 = squareInfo.integrations2;
integrations3 = squareInfo.integrations3;
surfaceArea1  = squareInfo.surfaceArea1;
surfaceArea2  = squareInfo.surfaceArea2;
surfaceArea3  = squareInfo.surfaceArea3;
numSubSq1     = squareInfo.numSubSq1;
numSubSq2     = squareInfo.numSubSq2;
numSubSq3     = squareInfo.numSubSq3;

%% Generate quadrature for face 1 (y-z plane)

% Initialize storage
xPos1      = cell(numSubSq1, 1);
yPos1      = cell(numSubSq1, 1);
gamma1     = cell(numSubSq1, 1);
theta1     = cell(numSubSq1, 1);
weights1   = cell(numSubSq1, 1);
totWeights = 0;
weightMat  = zeros(4, 4);

% Begin secant method
maxIter = 100;

% Go through each sub-square
for iSub = 1 : numSubSq1
    
    fprintf('sub-square %i of %i for face 1 \n', iSub, numSubSq1);
    
    % Initial ratio for current sub-square
    rhoInit = rhoInitAll1(iSub);
    
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
                xPosStore = [midSubX1(iSub) + (maxSubSubX1{iSub}(1) - minSubSubX1{iSub}(1)) * rhoPrev(1),...
                             midSubX1(iSub) + (maxSubSubX1{iSub}(2) - minSubSubX1{iSub}(2)) * rhoPrev(2),...
                             midSubX1(iSub) - (maxSubSubX1{iSub}(3) - minSubSubX1{iSub}(3)) * rhoPrev(3),...
                             midSubX1(iSub) - (maxSubSubX1{iSub}(4) - minSubSubX1{iSub}(4)) * rhoInit];
                yPosStore = [midSubY1(iSub) + (maxSubSubY1{iSub}(1) - minSubSubY1{iSub}(1)) * rhoPrev(1),...
                             midSubY1(iSub) - (maxSubSubY1{iSub}(2) - minSubSubY1{iSub}(2)) * rhoPrev(2),...
                             midSubY1(iSub) - (maxSubSubY1{iSub}(3) - minSubSubY1{iSub}(3)) * rhoPrev(3),...
                             midSubY1(iSub) + (maxSubSubY1{iSub}(4) - minSubSubY1{iSub}(4)) * rhoInit];
                [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
                constants = Basis(gammaStore, thetaStore);
                if constants == zeros(4)
                    convergedInner = 1;
                end
                weightsStore = Weights(integrations1{iSub}, constants);
                
                % Store weights in weight matrix
                for i = 1 : 4
                    weightMat(1, i) = weightsStore(i);
                end
                
                % Form secant matrix moving one ratio at a time
                numConv = 0;
                for m = 1 : 3
                    rhoTest = rhoPrev;
                    rhoTest(m) = rhoNext(m);
                    xPosTemp = [midSubX1(iSub) + (maxSubSubX1{iSub}(1) - minSubSubX1{iSub}(1)) * rhoTest(1),...
                                midSubX1(iSub) + (maxSubSubX1{iSub}(2) - minSubSubX1{iSub}(2)) * rhoTest(2),...
                                midSubX1(iSub) - (maxSubSubX1{iSub}(3) - minSubSubX1{iSub}(3)) * rhoTest(3),...
                                midSubX1(iSub) - (maxSubSubX1{iSub}(4) - minSubSubX1{iSub}(4)) * rhoInit];
                    yPosTemp = [midSubY1(iSub) + (maxSubSubY1{iSub}(1) - minSubSubY1{iSub}(1)) * rhoTest(1),...
                                midSubY1(iSub) - (maxSubSubY1{iSub}(2) - minSubSubY1{iSub}(2)) * rhoTest(2),...
                                midSubY1(iSub) - (maxSubSubY1{iSub}(3) - minSubSubY1{iSub}(3)) * rhoTest(3),...
                                midSubY1(iSub) + (maxSubSubY1{iSub}(4) - minSubSubY1{iSub}(4)) * rhoInit];
                    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                    constants = Basis(gammaTemp, thetaTemp);
                    if constants == zeros(4)
                        convergedInner = 1;
                    end
                    weightsTemp = Weights(integrations1{iSub}, constants);
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
                        difference(i) = surfaceArea1{iSub}(i) - weightMat(1, i);
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
                            difference(indexI) = surfaceArea1{iSub}(i) - weightMat(1, i);
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
                    difference = surfaceArea1{iSub}(unconvergedID) - weightMat(1, unconvergedID);
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
                    if abs(surfaceArea1{iSub}(1) - weightMat(1, 1)) / surfaceArea1{iSub}(1) < 1e-5 && ...
                       abs(surfaceArea1{iSub}(2) - weightMat(1, 2)) / surfaceArea1{iSub}(2) < 1e-5 && ...
                       abs(surfaceArea1{iSub}(3) - weightMat(1, 3)) / surfaceArea1{iSub}(3) < 1e-5 && ...
                       abs(surfaceArea1{iSub}(4) - weightMat(1, 4)) / surfaceArea1{iSub}(4) < 1e-5 && ...
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
    xPos1{iSub}    = xPosStore;
    yPos1{iSub}    = yPosStore;
    gamma1{iSub}   = gammaStore;
    theta1{iSub}   = thetaStore;
    weights1{iSub} = weightsStore;
    totWeights = totWeights + sum(weights1{iSub});
    
end

[min_val, min_pos] = min(factor_store);
fprintf('The minimum factor for face 1 is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(factor_store);
fprintf('The maximum factor for face 1 is: %E located in sub-square: %i \n', max_val, max_pos);

[min_val, min_pos] = min(tol_store);
fprintf('The minimum tol for face 1 is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(tol_store);
fprintf('The maximum tol for face 1 is: %E located in sub-square: %i \n', max_val, max_pos);

%% Generate quadrature for face 2 (x-z plane)

% Initialize storage
xPos2    = cell(numSubSq2, 1);
yPos2    = cell(numSubSq2, 1);
gamma2   = cell(numSubSq2, 1);
theta2   = cell(numSubSq2, 1);
weights2 = cell(numSubSq2, 1);
factor_store = [];
tol_store    = [];

% Begin secant method
maxIter = 100;

% Go through each sub-square
for iSub = 1 : numSubSq2
    
    fprintf('sub-square %i of %i for face 2 \n', iSub, numSubSq2);
    
    % Initial ratio for current sub-square
    rhoInit = rhoInitAll2(iSub);
    
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
                xPosStore = [midSubX2(iSub) + (maxSubSubX2{iSub}(1) - minSubSubX2{iSub}(1)) * rhoPrev(1),...
                             midSubX2(iSub) + (maxSubSubX2{iSub}(2) - minSubSubX2{iSub}(2)) * rhoPrev(2),...
                             midSubX2(iSub) - (maxSubSubX2{iSub}(3) - minSubSubX2{iSub}(3)) * rhoPrev(3),...
                             midSubX2(iSub) - (maxSubSubX2{iSub}(4) - minSubSubX2{iSub}(4)) * rhoInit];
                yPosStore = [midSubY2(iSub) + (maxSubSubY2{iSub}(1) - minSubSubY2{iSub}(1)) * rhoPrev(1),...
                             midSubY2(iSub) - (maxSubSubY2{iSub}(2) - minSubSubY2{iSub}(2)) * rhoPrev(2),...
                             midSubY2(iSub) - (maxSubSubY2{iSub}(3) - minSubSubY2{iSub}(3)) * rhoPrev(3),...
                             midSubY2(iSub) + (maxSubSubY2{iSub}(4) - minSubSubY2{iSub}(4)) * rhoInit];
                [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
                constants = Basis(gammaStore, thetaStore);
                if constants == zeros(4)
                    convergedInner = 1;
                end
                weightsStore = Weights(integrations2{iSub}, constants);
                
                % Store weights in weight matrix
                for i = 1 : 4
                    weightMat(1, i) = weightsStore(i);
                end
                
                % Form secant matrix moving one ratio at a time
                numConv = 0;
                for m = 1 : 3
                    rhoTest = rhoPrev;
                    rhoTest(m) = rhoNext(m);
                    xPosTemp = [midSubX2(iSub) + (maxSubSubX2{iSub}(1) - minSubSubX2{iSub}(1)) * rhoTest(1),...
                                midSubX2(iSub) + (maxSubSubX2{iSub}(2) - minSubSubX2{iSub}(2)) * rhoTest(2),...
                                midSubX2(iSub) - (maxSubSubX2{iSub}(3) - minSubSubX2{iSub}(3)) * rhoTest(3),...
                                midSubX2(iSub) - (maxSubSubX2{iSub}(4) - minSubSubX2{iSub}(4)) * rhoInit];
                    yPosTemp = [midSubY2(iSub) + (maxSubSubY2{iSub}(1) - minSubSubY2{iSub}(1)) * rhoTest(1),...
                                midSubY2(iSub) - (maxSubSubY2{iSub}(2) - minSubSubY2{iSub}(2)) * rhoTest(2),...
                                midSubY2(iSub) - (maxSubSubY2{iSub}(3) - minSubSubY2{iSub}(3)) * rhoTest(3),...
                                midSubY2(iSub) + (maxSubSubY2{iSub}(4) - minSubSubY2{iSub}(4)) * rhoInit];
                    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                    constants = Basis(gammaTemp, thetaTemp);
                    if constants == zeros(4)
                        convergedInner = 1;
                    end
                    weightsTemp = Weights(integrations2{iSub}, constants);
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
                        difference(i) = surfaceArea2{iSub}(i) - weightMat(1, i);
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
                            difference(indexI) = surfaceArea2{iSub}(i) - weightMat(1, i);
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
                    difference = surfaceArea2{iSub}(unconvergedID) - weightMat(1, unconvergedID);
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
                    if abs(surfaceArea2{iSub}(1) - weightMat(1, 1)) / surfaceArea2{iSub}(1) < 1e-5 && ...
                       abs(surfaceArea2{iSub}(2) - weightMat(1, 2)) / surfaceArea2{iSub}(2) < 1e-5 && ...
                       abs(surfaceArea2{iSub}(3) - weightMat(1, 3)) / surfaceArea2{iSub}(3) < 1e-5 && ...
                       abs(surfaceArea2{iSub}(4) - weightMat(1, 4)) / surfaceArea2{iSub}(4) < 1e-5 && ...
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
    xPos2{iSub}  = xPosStore;
    yPos2{iSub}  = yPosStore;
    gamma2{iSub} = gammaStore;
    theta2{iSub} = thetaStore;
    weights2{iSub} = weightsStore;
    totWeights = totWeights + sum(weights2{iSub});
    
end

[min_val, min_pos] = min(factor_store);
fprintf('The minimum factor for face 2 is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(factor_store);
fprintf('The maximum factor for face 2 is: %E located in sub-square: %i \n', max_val, max_pos);

[min_val, min_pos] = min(tol_store);
fprintf('The minimum tol for face 2 is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(tol_store);
fprintf('The maximum tol for face 2 is: %E located in sub-square: %i \n', max_val, max_pos);

% Rotate to face 2
for i = 1 : numSubSq2
    for j = 1 : 4        
        gamma2{i}(j) = (pi / 2) - atan(sqrt(3) * yPos2{i}(j));
        theta2{i}(j) = (pi / 2) - atan(sqrt(3) * xPos2{i}(j) / sqrt(1 + 3 * yPos2{i}(j) ^ 2));
    end
end

%% Generate quadrature for face 3 (x-y plane)

% Initialize storage
xPos3    = cell(numSubSq3, 1);
yPos3    = cell(numSubSq3, 1);
gamma3   = cell(numSubSq3, 1);
theta3   = cell(numSubSq3, 1);
weights3 = cell(numSubSq3, 1);
factor_store = [];
tol_store    = [];

% Begin secant method
maxIter = 100;

% Go through each sub-square
for iSub = 1 : numSubSq3
    
    fprintf('sub-square %i of %i for face 3 \n', iSub, numSubSq3);
    
    % Initial ratio for current sub-square
    rhoInit = rhoInitAll3(iSub);
    
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
                xPosStore = [midSubX3(iSub) + (maxSubSubX3{iSub}(1) - minSubSubX3{iSub}(1)) * rhoPrev(1),...
                             midSubX3(iSub) + (maxSubSubX3{iSub}(2) - minSubSubX3{iSub}(2)) * rhoPrev(2),...
                             midSubX3(iSub) - (maxSubSubX3{iSub}(3) - minSubSubX3{iSub}(3)) * rhoPrev(3),...
                             midSubX3(iSub) - (maxSubSubX3{iSub}(4) - minSubSubX3{iSub}(4)) * rhoInit];
                yPosStore = [midSubY3(iSub) + (maxSubSubY3{iSub}(1) - minSubSubY3{iSub}(1)) * rhoPrev(1),...
                             midSubY3(iSub) - (maxSubSubY3{iSub}(2) - minSubSubY3{iSub}(2)) * rhoPrev(2),...
                             midSubY3(iSub) - (maxSubSubY3{iSub}(3) - minSubSubY3{iSub}(3)) * rhoPrev(3),...
                             midSubY3(iSub) + (maxSubSubY3{iSub}(4) - minSubSubY3{iSub}(4)) * rhoInit];
                [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
                constants = Basis(gammaStore, thetaStore);
                if constants == zeros(4)
                    convergedInner = 1;
                end
                weightsStore = Weights(integrations3{iSub}, constants);
                
                % Store weights in weight matrix
                for i = 1 : 4
                    weightMat(1, i) = weightsStore(i);
                end
                
                % Form secant matrix moving one ratio at a time
                numConv = 0;
                for m = 1 : 3
                    rhoTest = rhoPrev;
                    rhoTest(m) = rhoNext(m);
                    xPosTemp = [midSubX3(iSub) + (maxSubSubX3{iSub}(1) - minSubSubX3{iSub}(1)) * rhoTest(1),...
                                midSubX3(iSub) + (maxSubSubX3{iSub}(2) - minSubSubX3{iSub}(2)) * rhoTest(2),...
                                midSubX3(iSub) - (maxSubSubX3{iSub}(3) - minSubSubX3{iSub}(3)) * rhoTest(3),...
                                midSubX3(iSub) - (maxSubSubX3{iSub}(4) - minSubSubX3{iSub}(4)) * rhoInit];
                    yPosTemp = [midSubY3(iSub) + (maxSubSubY3{iSub}(1) - minSubSubY3{iSub}(1)) * rhoTest(1),...
                                midSubY3(iSub) - (maxSubSubY3{iSub}(2) - minSubSubY3{iSub}(2)) * rhoTest(2),...
                                midSubY3(iSub) - (maxSubSubY3{iSub}(3) - minSubSubY3{iSub}(3)) * rhoTest(3),...
                                midSubY3(iSub) + (maxSubSubY3{iSub}(4) - minSubSubY3{iSub}(4)) * rhoInit];
                    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                    constants = Basis(gammaTemp, thetaTemp);
                    if constants == zeros(4)
                        convergedInner = 1;
                    end
                    weightsTemp = Weights(integrations3{iSub}, constants);
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
                        difference(i) = surfaceArea3{iSub}(i) - weightMat(1, i);
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
                            difference(indexI) = surfaceArea3{iSub}(i) - weightMat(1, i);
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
                    difference = surfaceArea3{iSub}(unconvergedID) - weightMat(1, unconvergedID);
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
                    if abs(surfaceArea3{iSub}(1) - weightMat(1, 1)) / surfaceArea3{iSub}(1) < 1e-5 && ...
                       abs(surfaceArea3{iSub}(2) - weightMat(1, 2)) / surfaceArea3{iSub}(2) < 1e-5 && ...
                       abs(surfaceArea3{iSub}(3) - weightMat(1, 3)) / surfaceArea3{iSub}(3) < 1e-5 && ...
                       abs(surfaceArea3{iSub}(4) - weightMat(1, 4)) / surfaceArea3{iSub}(4) < 1e-5 && ...
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
    xPos3{iSub}     = xPosStore;
    yPos3{iSub}     = yPosStore;
    weights3{iSub}  = weightsStore;
    totWeights = totWeights + sum(weights3{iSub});
    
end

[min_val, min_pos] = min(factor_store);
fprintf('The minimum factor for face 3 is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(factor_store);
fprintf('The maximum factor for face 3 is: %E located in sub-square: %i \n', max_val, max_pos);

[min_val, min_pos] = min(tol_store);
fprintf('The minimum tol for face 3 is: %E located in sub-square: %i \n', min_val, min_pos);
[max_val, max_pos] = max(tol_store);
fprintf('The maximum tol for face 3 is: %E located in sub-square: %i \n', max_val, max_pos);

% Rotate to face 3
for i = 1 : numSubSq3
    for j = 1 : 4        
        gamma3{i}(j) = atan(yPos3{i}(j) / xPos3{i}(j));
        theta3{i}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos3{i}(j) ^ 2 + xPos3{i}(j) ^ 2)));
    end
end

%% Maximum surface area error
test = zeros((numSubSq1 + numSubSq2 + numSubSq3) * 4, 1);
counter = 1;
for i = 1 : numSubSq1
    for j = 1 : 4      
        test(counter) = abs(weights1{i}(j) - surfaceArea1{i}(j)) / surfaceArea1{i}(j);
        counter = counter + 1;
    end
end
for i = 1 : numSubSq2
    for j = 1 : 4      
        test(counter) = abs(weights2{i}(j) - surfaceArea2{i}(j)) / surfaceArea2{i}(j);
        counter = counter + 1;
    end
end
for i = 1 : numSubSq3
    for j = 1 : 4      
        test(counter) = abs(weights3{i}(j) - surfaceArea3{i}(j)) / surfaceArea3{i}(j);
        counter = counter + 1;
    end
end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

%% Store LDFE-center quadrature data
quadrature.gamma1   = gamma1;
quadrature.gamma2   = gamma2;
quadrature.gamma3   = gamma3;
quadrature.theta1   = theta1;
quadrature.theta2   = theta2;
quadrature.theta3   = theta3;
quadrature.weights1 = weights1;
quadrature.weights2 = weights2;
quadrature.weights3 = weights3;

end