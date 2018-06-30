% Generate QDFE-ratio quadrature
function quadrature = QDFEsa(squareInfo)

%% Retrieve square properties
minSubSubX1    = squareInfo.minSubSubX1;
maxSubSubX1    = squareInfo.maxSubSubX1;
minSubSubY1    = squareInfo.minSubSubY1;
maxSubSubY1    = squareInfo.maxSubSubY1;
subSubLengthX1 = squareInfo.subSubLengthX1;
subSubLengthY1 = squareInfo.subSubLengthY1;
surfaceArea1   = squareInfo.surfaceArea1;
integrations1  = squareInfo.integrations1;
numSubSq1      = squareInfo.numSubSq1;

minSubSubX2    = squareInfo.minSubSubX2;
maxSubSubX2    = squareInfo.maxSubSubX2;
minSubSubY2    = squareInfo.minSubSubY2;
maxSubSubY2    = squareInfo.maxSubSubY2;
subSubLengthX2 = squareInfo.subSubLengthX2;
subSubLengthY2 = squareInfo.subSubLengthY2;
surfaceArea2   = squareInfo.surfaceArea2;
integrations2  = squareInfo.integrations2;
numSubSq2      = squareInfo.numSubSq2;

minSubSubX3    = squareInfo.minSubSubX3;
maxSubSubX3    = squareInfo.maxSubSubX3;
minSubSubY3    = squareInfo.minSubSubY3;
maxSubSubY3    = squareInfo.maxSubSubY3;
subSubLengthX3 = squareInfo.subSubLengthX3;
subSubLengthY3 = squareInfo.subSubLengthY3;
surfaceArea3   = squareInfo.surfaceArea3;
integrations3  = squareInfo.integrations3;
numSubSq3      = squareInfo.numSubSq3;

%% Initial ratio for QDFE-ratio method
[rho_corner_all1, rho_side_all1, ...
 rho_corner_all2, rho_side_all2, ...
 rho_corner_all3, rho_side_all3] = InitialRatio(squareInfo);

%% Generate quadrature for face 1 (y-z plane)

% Initialize storage
xPos1     = cell(numSubSq1, 1);
yPos1     = cell(numSubSq1, 1);
gamma1    = cell(numSubSq1, 1);
theta1    = cell(numSubSq1, 1);
weights1  = cell(numSubSq1, 1);
weightMat = zeros(8, 8);

% Begin secant method
maxIter = 100;

% Go through each sub-square
factorCounter = 1;
for iSub = 1 : numSubSq1
        
    fprintf('face1: sub-square %i of %i \n', iSub, numSubSq1);
    
    % Surface area without center sub-sub-square
    surfaceAreaTemp(1) = surfaceArea1{iSub}(1);
    surfaceAreaTemp(2) = surfaceArea1{iSub}(2);
    surfaceAreaTemp(3) = surfaceArea1{iSub}(3);
    surfaceAreaTemp(4) = surfaceArea1{iSub}(4);
    surfaceAreaTemp(5) = surfaceArea1{iSub}(6);
    surfaceAreaTemp(6) = surfaceArea1{iSub}(7);
    surfaceAreaTemp(7) = surfaceArea1{iSub}(8);
    surfaceAreaTemp(8) = surfaceArea1{iSub}(9);
    
    % Initial ratio for current sub-square
    rho_corner = rho_corner_all1(iSub);
    rho_side   = rho_side_all1(iSub);
    
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
                 maxSubSubX1{iSub}(1) - subSubLengthX1{iSub}(1) * rhoPrev(1), ...
                (maxSubSubX1{iSub}(2) + minSubSubX1{iSub}(2)) / 2, ...
                 minSubSubX1{iSub}(3) + subSubLengthX1{iSub}(3) * rhoPrev(3), ...
                 ...
                 maxSubSubX1{iSub}(4) - subSubLengthX1{iSub}(4) * rhoPrev(4), ...
                (maxSubSubX1{iSub}(5) + minSubSubX1{iSub}(5)) / 2, ...
                 minSubSubX1{iSub}(6) + subSubLengthX1{iSub}(6) * rhoPrev(5), ...
                 ...
                 maxSubSubX1{iSub}(7) - subSubLengthX1{iSub}(7) * rhoPrev(6), ...
                (maxSubSubX1{iSub}(8) + minSubSubX1{iSub}(8)) / 2, ...
                 minSubSubX1{iSub}(9) + subSubLengthX1{iSub}(9) * rhoPrev(8)];
             
            yPosStore = [...
                 maxSubSubY1{iSub}(1) - subSubLengthY1{iSub}(1) * rhoPrev(1), ...
                 maxSubSubY1{iSub}(2) - subSubLengthY1{iSub}(2) * rhoPrev(2), ...
                 maxSubSubY1{iSub}(3) - subSubLengthY1{iSub}(3) * rhoPrev(3), ...
                 ...
                (maxSubSubY1{iSub}(4) + minSubSubY1{iSub}(4)) / 2, ...
                (maxSubSubY1{iSub}(5) + minSubSubY1{iSub}(5)) / 2, ...
                (maxSubSubY1{iSub}(6) + minSubSubY1{iSub}(6)) / 2, ...
                ...
                 minSubSubY1{iSub}(7) + subSubLengthY1{iSub}(7) * rhoPrev(6), ...
                 minSubSubY1{iSub}(8) + subSubLengthY1{iSub}(8) * rhoPrev(7), ...
                 minSubSubY1{iSub}(9) + subSubLengthY1{iSub}(9) * rhoPrev(8)];
            [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
            constants = Basis(gammaStore, thetaStore);
            if constants == zeros(9)
                convergedInner = 1;
            end
            weightsStore = Weights(integrations1{iSub}, constants);
            
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
                     maxSubSubX1{iSub}(1) - subSubLengthX1{iSub}(1) * rhoTest(1), ...
                    (maxSubSubX1{iSub}(2) + minSubSubX1{iSub}(2)) / 2, ...
                     minSubSubX1{iSub}(3) + subSubLengthX1{iSub}(3) * rhoTest(3), ...
                     ...
                     maxSubSubX1{iSub}(4) - subSubLengthX1{iSub}(4) * rhoTest(4), ...
                    (maxSubSubX1{iSub}(5) + minSubSubX1{iSub}(5)) / 2, ...
                     minSubSubX1{iSub}(6) + subSubLengthX1{iSub}(6) * rhoTest(5), ...
                     ...
                     maxSubSubX1{iSub}(7) - subSubLengthX1{iSub}(7) * rhoTest(6), ...
                    (maxSubSubX1{iSub}(8) + minSubSubX1{iSub}(8)) / 2, ...
                     minSubSubX1{iSub}(9) + subSubLengthX1{iSub}(9) * rhoTest(8)];
                 
                yPosTemp = [...
                     maxSubSubY1{iSub}(1) - subSubLengthY1{iSub}(1) * rhoTest(1), ...
                     maxSubSubY1{iSub}(2) - subSubLengthY1{iSub}(2) * rhoTest(2), ...
                     maxSubSubY1{iSub}(3) - subSubLengthY1{iSub}(3) * rhoTest(3), ...
                     ...
                    (maxSubSubY1{iSub}(4) + minSubSubY1{iSub}(4)) / 2, ...
                    (maxSubSubY1{iSub}(5) + minSubSubY1{iSub}(5)) / 2, ...
                    (maxSubSubY1{iSub}(6) + minSubSubY1{iSub}(6)) / 2, ...
                    ...
                     minSubSubY1{iSub}(7) + subSubLengthY1{iSub}(7) * rhoTest(6), ...
                     minSubSubY1{iSub}(8) + subSubLengthY1{iSub}(8) * rhoTest(7), ...
                     minSubSubY1{iSub}(9) + subSubLengthY1{iSub}(9) * rhoTest(8)];
                [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                constants = Basis(gammaTemp, thetaTemp);
                if constants == zeros(9)
                    convergedInner = 1;
                end
                weightsTemp = Weights(integrations1{iSub}, constants);
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
    factor_store(factorCounter) = factor * 10;
    
    % Store error of each sub-square
    tol_store(factorCounter) = tol;
    
    % Check if current sub-square failed
    if convergedOuter == 0;
        error('QDFE-sa failed!');
    end
    
    % Store quadrature
    xPos1{iSub}     = xPosStore;
    yPos1{iSub}     = yPosStore;
    gamma1{iSub, 1} = gammaStore;
    theta1{iSub, 1} = thetaStore;
    weights1{iSub}  = weightsStore;
    
    % Update factor counter
    factorCounter = factorCounter + 1;
    
end

%% Generate quadrature for face 2 (x-z plane)

% Initialize storage
xPos2     = cell(numSubSq2, 1);
yPos2     = cell(numSubSq2, 1);
gamma2    = cell(numSubSq2, 1);
theta2    = cell(numSubSq2, 1);
weights2  = cell(numSubSq2, 1);
weightMat = zeros(8, 8);

% Begin secant method
maxIter = 100;

% Go through each sub-square
for iSub = 1 : numSubSq2
        
    fprintf('face2: sub-square %i of %i \n', iSub, numSubSq2);

    % Surface area without center sub-sub-square
    surfaceAreaTemp(1) = surfaceArea2{iSub}(1);
    surfaceAreaTemp(2) = surfaceArea2{iSub}(2);
    surfaceAreaTemp(3) = surfaceArea2{iSub}(3);
    surfaceAreaTemp(4) = surfaceArea2{iSub}(4);
    surfaceAreaTemp(5) = surfaceArea2{iSub}(6);
    surfaceAreaTemp(6) = surfaceArea2{iSub}(7);
    surfaceAreaTemp(7) = surfaceArea2{iSub}(8);
    surfaceAreaTemp(8) = surfaceArea2{iSub}(9);
    
    % Initial ratio for current sub-square
    rho_corner = rho_corner_all2(iSub);
    rho_side   = rho_side_all2(iSub);
    
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
                 maxSubSubX2{iSub}(1) - subSubLengthX2{iSub}(1) * rhoPrev(1), ...
                (maxSubSubX2{iSub}(2) + minSubSubX2{iSub}(2)) / 2, ...
                 minSubSubX2{iSub}(3) + subSubLengthX2{iSub}(3) * rhoPrev(3), ...
                 ...
                 maxSubSubX2{iSub}(4) - subSubLengthX2{iSub}(4) * rhoPrev(4), ...
                (maxSubSubX2{iSub}(5) + minSubSubX2{iSub}(5)) / 2, ...
                 minSubSubX2{iSub}(6) + subSubLengthX2{iSub}(6) * rhoPrev(5), ...
                 ...
                 maxSubSubX2{iSub}(7) - subSubLengthX2{iSub}(7) * rhoPrev(6), ...
                (maxSubSubX2{iSub}(8) + minSubSubX2{iSub}(8)) / 2, ...
                 minSubSubX2{iSub}(9) + subSubLengthX2{iSub}(9) * rhoPrev(8)];
             
            yPosStore = [...
                 maxSubSubY2{iSub}(1) - subSubLengthY2{iSub}(1) * rhoPrev(1), ...
                 maxSubSubY2{iSub}(2) - subSubLengthY2{iSub}(2) * rhoPrev(2), ...
                 maxSubSubY2{iSub}(3) - subSubLengthY2{iSub}(3) * rhoPrev(3), ...
                 ...
                (maxSubSubY2{iSub}(4) + minSubSubY2{iSub}(4)) / 2, ...
                (maxSubSubY2{iSub}(5) + minSubSubY2{iSub}(5)) / 2, ...
                (maxSubSubY2{iSub}(6) + minSubSubY2{iSub}(6)) / 2, ...
                ...
                 minSubSubY2{iSub}(7) + subSubLengthY2{iSub}(7) * rhoPrev(6), ...
                 minSubSubY2{iSub}(8) + subSubLengthY2{iSub}(8) * rhoPrev(7), ...
                 minSubSubY2{iSub}(9) + subSubLengthY2{iSub}(9) * rhoPrev(8)];
            [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
            constants = Basis(gammaStore, thetaStore);
            if constants == zeros(9)
                convergedInner = 1;
            end
            weightsStore = Weights(integrations2{iSub}, constants);
            
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
                     maxSubSubX2{iSub}(1) - subSubLengthX2{iSub}(1) * rhoTest(1), ...
                    (maxSubSubX2{iSub}(2) + minSubSubX2{iSub}(2)) / 2, ...
                     minSubSubX2{iSub}(3) + subSubLengthX2{iSub}(3) * rhoTest(3), ...
                     ...
                     maxSubSubX2{iSub}(4) - subSubLengthX2{iSub}(4) * rhoTest(4), ...
                    (maxSubSubX2{iSub}(5) + minSubSubX2{iSub}(5)) / 2, ...
                     minSubSubX2{iSub}(6) + subSubLengthX2{iSub}(6) * rhoTest(5), ...
                     ...
                     maxSubSubX2{iSub}(7) - subSubLengthX2{iSub}(7) * rhoTest(6), ...
                    (maxSubSubX2{iSub}(8) + minSubSubX2{iSub}(8)) / 2, ...
                     minSubSubX2{iSub}(9) + subSubLengthX2{iSub}(9) * rhoTest(8)];
                 
                yPosTemp = [...
                     maxSubSubY2{iSub}(1) - subSubLengthY2{iSub}(1) * rhoTest(1), ...
                     maxSubSubY2{iSub}(2) - subSubLengthY2{iSub}(2) * rhoTest(2), ...
                     maxSubSubY2{iSub}(3) - subSubLengthY2{iSub}(3) * rhoTest(3), ...
                     ...
                    (maxSubSubY2{iSub}(4) + minSubSubY2{iSub}(4)) / 2, ...
                    (maxSubSubY2{iSub}(5) + minSubSubY2{iSub}(5)) / 2, ...
                    (maxSubSubY2{iSub}(6) + minSubSubY2{iSub}(6)) / 2, ...
                    ...
                     minSubSubY2{iSub}(7) + subSubLengthY2{iSub}(7) * rhoTest(6), ...
                     minSubSubY2{iSub}(8) + subSubLengthY2{iSub}(8) * rhoTest(7), ...
                     minSubSubY2{iSub}(9) + subSubLengthY2{iSub}(9) * rhoTest(8)];
                [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                constants = Basis(gammaTemp, thetaTemp);
                if constants == zeros(9)
                    convergedInner = 1;
                end
                weightsTemp = Weights(integrations2{iSub}, constants);
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
    factor_store(factorCounter) = factor * 10;
    
    % Store error of each sub-square
    tol_store(factorCounter) = tol;
    
    % Check if current sub-square failed
    if convergedOuter == 0;
        error('QDFE-sa failed!');
    end
    
    % Store quadrature
    xPos2{iSub}     = xPosStore;
    yPos2{iSub}     = yPosStore;
    gamma2{iSub, 1} = gammaStore;
    theta2{iSub, 1} = thetaStore;
    weights2{iSub}  = weightsStore;
    
    % Update factor counter
    factorCounter = factorCounter + 1;
    
end

% Rotate to face 2
for i = 1 : numSubSq2
    for j = 1 : 9        
        gamma2{i}(j) = (pi / 2) - atan(sqrt(3) * yPos2{i}(j));
        theta2{i}(j) = (pi / 2) - atan(sqrt(3) * xPos2{i}(j) / sqrt(1 + 3 * yPos2{i}(j) ^ 2));
    end
end

%% Generate quadrature for face 3 (x-y plane)

% Initialize storage
xPos3     = cell(numSubSq3, 1);
yPos3     = cell(numSubSq3, 1);
gamma3    = cell(numSubSq3, 1);
theta3    = cell(numSubSq3, 1);
weights3  = cell(numSubSq3, 1);
weightMat = zeros(8, 8);

% Begin secant method
maxIter = 100;

% Go through each sub-square
for iSub = 1 : numSubSq3
        
    fprintf('face3: sub-square %i of %i \n', iSub, numSubSq3);
        
    % Surface area without center sub-sub-square
    surfaceAreaTemp(1) = surfaceArea3{iSub}(1);
    surfaceAreaTemp(2) = surfaceArea3{iSub}(2);
    surfaceAreaTemp(3) = surfaceArea3{iSub}(3);
    surfaceAreaTemp(4) = surfaceArea3{iSub}(4);
    surfaceAreaTemp(5) = surfaceArea3{iSub}(6);
    surfaceAreaTemp(6) = surfaceArea3{iSub}(7);
    surfaceAreaTemp(7) = surfaceArea3{iSub}(8);
    surfaceAreaTemp(8) = surfaceArea3{iSub}(9);
    
    % Initial ratio for current sub-square
    rho_corner = rho_corner_all3(iSub);
    rho_side   = rho_side_all3(iSub);
    
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
                 maxSubSubX3{iSub}(1) - subSubLengthX3{iSub}(1) * rhoPrev(1), ...
                (maxSubSubX3{iSub}(2) + minSubSubX3{iSub}(2)) / 2, ...
                 minSubSubX3{iSub}(3) + subSubLengthX3{iSub}(3) * rhoPrev(3), ...
                 ...
                 maxSubSubX3{iSub}(4) - subSubLengthX3{iSub}(4) * rhoPrev(4), ...
                (maxSubSubX3{iSub}(5) + minSubSubX3{iSub}(5)) / 2, ...
                 minSubSubX3{iSub}(6) + subSubLengthX3{iSub}(6) * rhoPrev(5), ...
                 ...
                 maxSubSubX3{iSub}(7) - subSubLengthX3{iSub}(7) * rhoPrev(6), ...
                (maxSubSubX3{iSub}(8) + minSubSubX3{iSub}(8)) / 2, ...
                 minSubSubX3{iSub}(9) + subSubLengthX3{iSub}(9) * rhoPrev(8)];
             
            yPosStore = [...
                 maxSubSubY3{iSub}(1) - subSubLengthY3{iSub}(1) * rhoPrev(1), ...
                 maxSubSubY3{iSub}(2) - subSubLengthY3{iSub}(2) * rhoPrev(2), ...
                 maxSubSubY3{iSub}(3) - subSubLengthY3{iSub}(3) * rhoPrev(3), ...
                 ...
                (maxSubSubY3{iSub}(4) + minSubSubY3{iSub}(4)) / 2, ...
                (maxSubSubY3{iSub}(5) + minSubSubY3{iSub}(5)) / 2, ...
                (maxSubSubY3{iSub}(6) + minSubSubY3{iSub}(6)) / 2, ...
                ...
                 minSubSubY3{iSub}(7) + subSubLengthY3{iSub}(7) * rhoPrev(6), ...
                 minSubSubY3{iSub}(8) + subSubLengthY3{iSub}(8) * rhoPrev(7), ...
                 minSubSubY3{iSub}(9) + subSubLengthY3{iSub}(9) * rhoPrev(8)];
            [gammaStore, thetaStore] = LocalToGlobal(xPosStore, yPosStore);
            constants = Basis(gammaStore, thetaStore);
            if constants == zeros(9)
                convergedInner = 1;
            end
            weightsStore = Weights(integrations3{iSub}, constants);
            
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
                     maxSubSubX3{iSub}(1) - subSubLengthX3{iSub}(1) * rhoTest(1), ...
                    (maxSubSubX3{iSub}(2) + minSubSubX3{iSub}(2)) / 2, ...
                     minSubSubX3{iSub}(3) + subSubLengthX3{iSub}(3) * rhoTest(3), ...
                     ...
                     maxSubSubX3{iSub}(4) - subSubLengthX3{iSub}(4) * rhoTest(4), ...
                    (maxSubSubX3{iSub}(5) + minSubSubX3{iSub}(5)) / 2, ...
                     minSubSubX3{iSub}(6) + subSubLengthX3{iSub}(6) * rhoTest(5), ...
                     ...
                     maxSubSubX3{iSub}(7) - subSubLengthX3{iSub}(7) * rhoTest(6), ...
                    (maxSubSubX3{iSub}(8) + minSubSubX3{iSub}(8)) / 2, ...
                     minSubSubX3{iSub}(9) + subSubLengthX3{iSub}(9) * rhoTest(8)];
                 
                yPosTemp = [...
                     maxSubSubY3{iSub}(1) - subSubLengthY3{iSub}(1) * rhoTest(1), ...
                     maxSubSubY3{iSub}(2) - subSubLengthY3{iSub}(2) * rhoTest(2), ...
                     maxSubSubY3{iSub}(3) - subSubLengthY3{iSub}(3) * rhoTest(3), ...
                     ...
                    (maxSubSubY3{iSub}(4) + minSubSubY3{iSub}(4)) / 2, ...
                    (maxSubSubY3{iSub}(5) + minSubSubY3{iSub}(5)) / 2, ...
                    (maxSubSubY3{iSub}(6) + minSubSubY3{iSub}(6)) / 2, ...
                    ...
                     minSubSubY3{iSub}(7) + subSubLengthY3{iSub}(7) * rhoTest(6), ...
                     minSubSubY3{iSub}(8) + subSubLengthY3{iSub}(8) * rhoTest(7), ...
                     minSubSubY3{iSub}(9) + subSubLengthY3{iSub}(9) * rhoTest(8)];
                [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
                constants = Basis(gammaTemp, thetaTemp);
                if constants == zeros(9)
                    convergedInner = 1;
                end
                weightsTemp = Weights(integrations3{iSub}, constants);
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
    factor_store(factorCounter) = factor * 10;
    
    % Store error of each sub-square
    tol_store(factorCounter) = tol;
    
    % Check if current sub-square failed
    if convergedOuter == 0;
        error('QDFE-sa failed!');
    end
    
    % Store quadrature
    xPos3{iSub}     = xPosStore;
    yPos3{iSub}     = yPosStore;
    gamma3{iSub, 1} = gammaStore;
    theta3{iSub, 1} = thetaStore;
    weights3{iSub}  = weightsStore;
    
    % Update factor counter
    factorCounter = factorCounter + 1;
    
end

% Rotate to face 3
for i = 1 : numSubSq3
    for j = 1 : 9
        gamma3{i}(j) = atan(yPos3{i}(j) / xPos3{i}(j));
        theta3{i}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos3{i}(j) ^ 2 + xPos3{i}(j) ^ 2)));
    end
end

%% Maximum surface area error
test = zeros((numSubSq1 + numSubSq2 + numSubSq3) * 9, 1);
counter = 1;
for i = 1 : numSubSq1
    for j = 1 : 9
        test(counter) = abs(weights1{i}(j) - surfaceArea1{i}(j)) / surfaceArea1{i}(j);
        counter = counter + 1;        
    end

end
for i = 1 : numSubSq2
    for j = 1 : 9
        test(counter) = abs(weights2{i}(j) - surfaceArea2{i}(j)) / surfaceArea2{i}(j);
        counter = counter + 1;        
    end

end
for i = 1 : numSubSq3
    for j = 1 : 9
        test(counter) = abs(weights3{i}(j) - surfaceArea3{i}(j)) / surfaceArea3{i}(j);
        counter = counter + 1;        
    end

end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

%% Store QDFE-sa quadrature data
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