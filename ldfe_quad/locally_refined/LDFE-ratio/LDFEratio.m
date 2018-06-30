% Generate LDFE-ratio quadrature
function quadrature = LDFEratio(squareInfo, id)

%% Retrieve square properties
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
xPos1    = cell(numSubSq1, 1);
yPos1    = cell(numSubSq1, 1);
gamma1   = cell(numSubSq1, 1);
theta1   = cell(numSubSq1, 1);
weights1 = cell(numSubSq1, 1);
totWeights = 0;

% Define iteration parameters
areaEps   = 1e-15;
maxCounts = 10000;

% Go through each sub-square
for iSub = 1 : numSubSq1
    
    % Upper bound guess
    ratio_a = 0.8;
    
    % Lower bound guess
    ratio_b = 0.4;
    
    % Solve initial weights for a ratio
    xPosTemp = [...
        midSubX1(iSub) + (maxSubSubX1{iSub}(1) - minSubSubX1{iSub}(1)) * ratio_a,...
        midSubX1(iSub) + (maxSubSubX1{iSub}(2) - minSubSubX1{iSub}(2)) * ratio_a,...
        midSubX1(iSub) - (maxSubSubX1{iSub}(3) - minSubSubX1{iSub}(3)) * ratio_a,...
        midSubX1(iSub) - (maxSubSubX1{iSub}(4) - minSubSubX1{iSub}(4)) * ratio_a];
    yPosTemp = [...
        midSubY1(iSub) + (maxSubSubY1{iSub}(1) - minSubSubY1{iSub}(1)) * ratio_a,...
        midSubY1(iSub) - (maxSubSubY1{iSub}(2) - minSubSubY1{iSub}(2)) * ratio_a,...
        midSubY1(iSub) - (maxSubSubY1{iSub}(3) - minSubSubY1{iSub}(3)) * ratio_a,...
        midSubY1(iSub) + (maxSubSubY1{iSub}(4) - minSubSubY1{iSub}(4)) * ratio_a];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations1{iSub}, constants);
    
    % Initial residual
    res_a = (surfaceArea1{iSub}(id) - weightsTemp(id)) / surfaceArea1{iSub}(id);
    
    % Make sure residual is positive
    if res_a < 0
        error('ratio_a gave negative weights!');
    end
    
    % Solve new weights
    xPosTemp = [...
        midSubX1(iSub) + (maxSubSubX1{iSub}(1) - minSubSubX1{iSub}(1)) * ratio_b,...
        midSubX1(iSub) + (maxSubSubX1{iSub}(2) - minSubSubX1{iSub}(2)) * ratio_b,...
        midSubX1(iSub) - (maxSubSubX1{iSub}(3) - minSubSubX1{iSub}(3)) * ratio_b,...
        midSubX1(iSub) - (maxSubSubX1{iSub}(4) - minSubSubX1{iSub}(4)) * ratio_b];
    yPosTemp = [...
        midSubY1(iSub) + (maxSubSubY1{iSub}(1) - minSubSubY1{iSub}(1)) * ratio_b,...
        midSubY1(iSub) - (maxSubSubY1{iSub}(2) - minSubSubY1{iSub}(2)) * ratio_b,...
        midSubY1(iSub) - (maxSubSubY1{iSub}(3) - minSubSubY1{iSub}(3)) * ratio_b,...
        midSubY1(iSub) + (maxSubSubY1{iSub}(4) - minSubSubY1{iSub}(4)) * ratio_b];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations1{iSub}, constants);
    
    % Solve new residual
    res_b = (surfaceArea1{iSub}(id) - weightsTemp(id)) / surfaceArea1{iSub}(id);
    
    % Make sure residual is negative
    if res_b > 0
        error('ratio_b is positive!');
    end
    
    % Iterate until chosen sub-square weight equals its surface area
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts                
        
        % Calculate next guess
        x = (ratio_b + ratio_a) / 2;
        
        % New weights for next guess
        xPosTemp = [midSubX1(iSub) + (maxSubSubX1{iSub}(1) - minSubSubX1{iSub}(1)) * x,...
                    midSubX1(iSub) + (maxSubSubX1{iSub}(2) - minSubSubX1{iSub}(2)) * x,...
                    midSubX1(iSub) - (maxSubSubX1{iSub}(3) - minSubSubX1{iSub}(3)) * x,...
                    midSubX1(iSub) - (maxSubSubX1{iSub}(4) - minSubSubX1{iSub}(4)) * x];
        yPosTemp = [midSubY1(iSub) + (maxSubSubY1{iSub}(1) - minSubSubY1{iSub}(1)) * x,...
                    midSubY1(iSub) - (maxSubSubY1{iSub}(2) - minSubSubY1{iSub}(2)) * x,...
                    midSubY1(iSub) - (maxSubSubY1{iSub}(3) - minSubSubY1{iSub}(3)) * x,...
                    midSubY1(iSub) + (maxSubSubY1{iSub}(4) - minSubSubY1{iSub}(4)) * x];
        [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
        constants = Basis(gammaTemp, thetaTemp);
        weightsTemp = Weights(integrations1{iSub}, constants);
        
        % New residual
        res_x = (surfaceArea1{iSub}(id) - weightsTemp(id)) / surfaceArea1{iSub}(id);
                
        % Ensure residuals are not equal
        if res_a == res_b
            error('Error: residuals are equal');
        else            
            
            % Check for convergence
            if abs(ratio_a - ratio_b) < areaEps               
                converged = 1;            
            else
                
                % If the residual of next guess has same sign as a
                if res_x / res_a > 0
                    
                    % New a ratio
                    ratio_a = x;
                    
                    % New a residual
                    res_a = res_x;
                    
                else
                    % New b ratio
                    ratio_b = x;
                    
                    % New b residual
                    res_b = res_x;
                    
                end
                
                % Update counter
                counter = counter + 1;
                
            end
        end
    end
    
    % Check for failure
    if converged == 0
        error('Error: LDFE-ratio failed');
    end
    
    % Store quadrature data
    xPos1{iSub}    = xPosTemp;
    yPos1{iSub}    = yPosTemp;
    gamma1{iSub}   = gammaTemp;
    theta1{iSub}   = thetaTemp;
    weights1{iSub} = weightsTemp;
    totWeights = totWeights + sum(weights1{iSub});
end

%% Generate quadrature for face 2 (x-z plane)

% Initialize storage
xPos2    = cell(numSubSq2, 1);
yPos2    = cell(numSubSq2, 1);
gamma2   = cell(numSubSq2, 1);
theta2   = cell(numSubSq2, 1);
weights2 = cell(numSubSq2, 1);

% Go through each sub-square
for iSub = 1 : numSubSq2
    
    % Upper bound guess
    ratio_a = 0.8;
    
    % Lower bound guess
    ratio_b = 0.4;
    
    % Solve initial weights for a ratio
    xPosTemp = [...
        midSubX2(iSub) + (maxSubSubX2{iSub}(1) - minSubSubX2{iSub}(1)) * ratio_a,...
        midSubX2(iSub) + (maxSubSubX2{iSub}(2) - minSubSubX2{iSub}(2)) * ratio_a,...
        midSubX2(iSub) - (maxSubSubX2{iSub}(3) - minSubSubX2{iSub}(3)) * ratio_a,...
        midSubX2(iSub) - (maxSubSubX2{iSub}(4) - minSubSubX2{iSub}(4)) * ratio_a];
    yPosTemp = [...
        midSubY2(iSub) + (maxSubSubY2{iSub}(1) - minSubSubY2{iSub}(1)) * ratio_a,...
        midSubY2(iSub) - (maxSubSubY2{iSub}(2) - minSubSubY2{iSub}(2)) * ratio_a,...
        midSubY2(iSub) - (maxSubSubY2{iSub}(3) - minSubSubY2{iSub}(3)) * ratio_a,...
        midSubY2(iSub) + (maxSubSubY2{iSub}(4) - minSubSubY2{iSub}(4)) * ratio_a];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations2{iSub}, constants);
    
    % Initial residual
    res_a = (surfaceArea2{iSub}(id) - weightsTemp(id)) / surfaceArea2{iSub}(id);
    
    % Make sure residual is positive
    if res_a < 0
        error('ratio_a gave negative weights!');
    end
    
    % Solve new weights
    xPosTemp = [...
        midSubX2(iSub) + (maxSubSubX2{iSub}(1) - minSubSubX2{iSub}(1)) * ratio_b,...
        midSubX2(iSub) + (maxSubSubX2{iSub}(2) - minSubSubX2{iSub}(2)) * ratio_b,...
        midSubX2(iSub) - (maxSubSubX2{iSub}(3) - minSubSubX2{iSub}(3)) * ratio_b,...
        midSubX2(iSub) - (maxSubSubX2{iSub}(4) - minSubSubX2{iSub}(4)) * ratio_b];
    yPosTemp = [...
        midSubY2(iSub) + (maxSubSubY2{iSub}(1) - minSubSubY2{iSub}(1)) * ratio_b,...
        midSubY2(iSub) - (maxSubSubY2{iSub}(2) - minSubSubY2{iSub}(2)) * ratio_b,...
        midSubY2(iSub) - (maxSubSubY2{iSub}(3) - minSubSubY2{iSub}(3)) * ratio_b,...
        midSubY2(iSub) + (maxSubSubY2{iSub}(4) - minSubSubY2{iSub}(4)) * ratio_b];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations2{iSub}, constants);
    
    % Solve new residual
    res_b = (surfaceArea2{iSub}(id) - weightsTemp(id)) / surfaceArea2{iSub}(id);
    
    % Make sure residual is negative
    if res_b > 0
        error('ratio_b is positive!');
    end
    
    % Iterate until chosen sub-square weight equals its surface area
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts                
        
        % Calculate next guess
        x = (ratio_b + ratio_a) / 2;
        
        % New weights for next guess
        xPosTemp = [midSubX2(iSub) + (maxSubSubX2{iSub}(1) - minSubSubX2{iSub}(1)) * x,...
                    midSubX2(iSub) + (maxSubSubX2{iSub}(2) - minSubSubX2{iSub}(2)) * x,...
                    midSubX2(iSub) - (maxSubSubX2{iSub}(3) - minSubSubX2{iSub}(3)) * x,...
                    midSubX2(iSub) - (maxSubSubX2{iSub}(4) - minSubSubX2{iSub}(4)) * x];
        yPosTemp = [midSubY2(iSub) + (maxSubSubY2{iSub}(1) - minSubSubY2{iSub}(1)) * x,...
                    midSubY2(iSub) - (maxSubSubY2{iSub}(2) - minSubSubY2{iSub}(2)) * x,...
                    midSubY2(iSub) - (maxSubSubY2{iSub}(3) - minSubSubY2{iSub}(3)) * x,...
                    midSubY2(iSub) + (maxSubSubY2{iSub}(4) - minSubSubY2{iSub}(4)) * x];
        [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
        constants = Basis(gammaTemp, thetaTemp);
        weightsTemp = Weights(integrations2{iSub}, constants);
        
        % New residual
        res_x = (surfaceArea2{iSub}(id) - weightsTemp(id)) / surfaceArea2{iSub}(id);
                
        % Ensure residuals are not equal
        if res_a == res_b
            error('Error: residuals are equal');
        else            
            
            % Check for convergence
            if abs(ratio_a - ratio_b) < areaEps               
                converged = 1;            
            else
                
                % If the residual of next guess has same sign as a
                if res_x / res_a > 0
                    
                    % New a ratio
                    ratio_a = x;
                    
                    % New a residual
                    res_a = res_x;
                    
                else
                    % New b ratio
                    ratio_b = x;
                    
                    % New b residual
                    res_b = res_x;
                    
                end
                
                % Update counter
                counter = counter + 1;
                
            end
        end
    end
    
    % Check for failure
    if converged == 0
        error('Error: LDFE-ratio failed');
    end
    
    % Store quadrature data
    xPos2{iSub}    = xPosTemp;
    yPos2{iSub}    = yPosTemp;        
    weights2{iSub} = weightsTemp;
    totWeights = totWeights + sum(weights2{iSub});
end

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

% Go through each sub-square
for iSub = 1 : numSubSq3
    
    % Upper bound guess
    ratio_a = 0.8;
    
    % Lower bound guess
    ratio_b = 0.4;
    
    % Solve initial weights for a ratio
    xPosTemp = [...
        midSubX3(iSub) + (maxSubSubX3{iSub}(1) - minSubSubX3{iSub}(1)) * ratio_a,...
        midSubX3(iSub) + (maxSubSubX3{iSub}(2) - minSubSubX3{iSub}(2)) * ratio_a,...
        midSubX3(iSub) - (maxSubSubX3{iSub}(3) - minSubSubX3{iSub}(3)) * ratio_a,...
        midSubX3(iSub) - (maxSubSubX3{iSub}(4) - minSubSubX3{iSub}(4)) * ratio_a];
    yPosTemp = [...
        midSubY3(iSub) + (maxSubSubY3{iSub}(1) - minSubSubY3{iSub}(1)) * ratio_a,...
        midSubY3(iSub) - (maxSubSubY3{iSub}(2) - minSubSubY3{iSub}(2)) * ratio_a,...
        midSubY3(iSub) - (maxSubSubY3{iSub}(3) - minSubSubY3{iSub}(3)) * ratio_a,...
        midSubY3(iSub) + (maxSubSubY3{iSub}(4) - minSubSubY3{iSub}(4)) * ratio_a];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations3{iSub}, constants);
    
    % Initial residual
    res_a = (surfaceArea3{iSub}(id) - weightsTemp(id)) / surfaceArea3{iSub}(id);
    
    % Make sure residual is positive
    if res_a < 0
        error('ratio_a gave negative weights!');
    end
    
    % Solve new weights
    xPosTemp = [...
        midSubX3(iSub) + (maxSubSubX3{iSub}(1) - minSubSubX3{iSub}(1)) * ratio_b,...
        midSubX3(iSub) + (maxSubSubX3{iSub}(2) - minSubSubX3{iSub}(2)) * ratio_b,...
        midSubX3(iSub) - (maxSubSubX3{iSub}(3) - minSubSubX3{iSub}(3)) * ratio_b,...
        midSubX3(iSub) - (maxSubSubX3{iSub}(4) - minSubSubX3{iSub}(4)) * ratio_b];
    yPosTemp = [...
        midSubY3(iSub) + (maxSubSubY3{iSub}(1) - minSubSubY3{iSub}(1)) * ratio_b,...
        midSubY3(iSub) - (maxSubSubY3{iSub}(2) - minSubSubY3{iSub}(2)) * ratio_b,...
        midSubY3(iSub) - (maxSubSubY3{iSub}(3) - minSubSubY3{iSub}(3)) * ratio_b,...
        midSubY3(iSub) + (maxSubSubY3{iSub}(4) - minSubSubY3{iSub}(4)) * ratio_b];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations3{iSub}, constants);
    
    % Solve new residual
    res_b = (surfaceArea3{iSub}(id) - weightsTemp(id)) / surfaceArea3{iSub}(id);
    
    % Make sure residual is negative
    if res_b > 0
        error('ratio_b is positive!');
    end
    
    % Iterate until chosen sub-square weight equals its surface area
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts                
        
        % Calculate next guess
        x = (ratio_b + ratio_a) / 2;
        
        % New weights for next guess
        xPosTemp = [midSubX3(iSub) + (maxSubSubX3{iSub}(1) - minSubSubX3{iSub}(1)) * x,...
                    midSubX3(iSub) + (maxSubSubX3{iSub}(2) - minSubSubX3{iSub}(2)) * x,...
                    midSubX3(iSub) - (maxSubSubX3{iSub}(3) - minSubSubX3{iSub}(3)) * x,...
                    midSubX3(iSub) - (maxSubSubX3{iSub}(4) - minSubSubX3{iSub}(4)) * x];
        yPosTemp = [midSubY3(iSub) + (maxSubSubY3{iSub}(1) - minSubSubY3{iSub}(1)) * x,...
                    midSubY3(iSub) - (maxSubSubY3{iSub}(2) - minSubSubY3{iSub}(2)) * x,...
                    midSubY3(iSub) - (maxSubSubY3{iSub}(3) - minSubSubY3{iSub}(3)) * x,...
                    midSubY3(iSub) + (maxSubSubY3{iSub}(4) - minSubSubY3{iSub}(4)) * x];
        [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
        constants = Basis(gammaTemp, thetaTemp);
        weightsTemp = Weights(integrations3{iSub}, constants);
        
        % New residual
        res_x = (surfaceArea3{iSub}(id) - weightsTemp(id)) / surfaceArea3{iSub}(id);
                
        % Ensure residuals are not equal
        if res_a == res_b
            error('Error: residuals are equal');
        else            
            
            % Check for convergence
            if abs(ratio_a - ratio_b) < areaEps               
                converged = 1;            
            else
                
                % If the residual of next guess has same sign as a
                if res_x / res_a > 0
                    
                    % New a ratio
                    ratio_a = x;
                    
                    % New a residual
                    res_a = res_x;
                    
                else
                    % New b ratio
                    ratio_b = x;
                    
                    % New b residual
                    res_b = res_x;
                    
                end
                
                % Update counter
                counter = counter + 1;
                
            end
        end
    end
    
    % Check for failure
    if converged == 0
        error('Error: LDFE-ratio failed');
    end
    
    % Store quadrature data
    xPos3{iSub}    = xPosTemp;
    yPos3{iSub}    = yPosTemp;
    weights3{iSub} = weightsTemp;
    totWeights = totWeights + sum(weights3{iSub});
end

% Rotate to face 3
for i = 1 : numSubSq3
    for j = 1 : 4        
        gamma3{i}(j) = atan(yPos3{i}(j) / xPos3{i}(j));
        theta3{i}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos3{i}(j) ^ 2 + xPos3{i}(j) ^ 2)));
    end
end

%% Maximum surface area error
test = zeros(numSubSq1 + numSubSq2 + numSubSq3, 1);
counter = 1;
for i = 1 : numSubSq1
    test(counter) = abs(weights1{i}(id) - surfaceArea1{i}(id)) / surfaceArea1{i}(id);
    counter = counter + 1;
end
for i = 1 : numSubSq2
    test(counter) = abs(weights2{i}(id) - surfaceArea2{i}(id)) / surfaceArea2{i}(id);
    counter = counter + 1;
end
for i = 1 : numSubSq3
    test(counter) = abs(weights3{i}(id) - surfaceArea3{i}(id)) / surfaceArea3{i}(id);
    counter = counter + 1;
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