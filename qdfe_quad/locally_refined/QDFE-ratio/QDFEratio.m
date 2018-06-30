% Generate QDFE-ratio quadrature
function quadrature = QDFEratio(squareInfo)

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
        
    % Find ratio_a (equal corner and side ratios) minimizing residual
    [ratio_0, res_a] = find_ratio_eq( ...
        subSubLengthX1{iSub}, subSubLengthY1{iSub}, ...
           minSubSubX1{iSub},    maxSubSubX1{iSub}, ...
           minSubSubY1{iSub},    maxSubSubY1{iSub}, ...
         integrations1{iSub},   surfaceArea1{iSub});
    ratio_a = [ratio_0, ratio_0];
    
    % Find ratio_b and ratio_c having opposite residual signs
    [ratio_corner, ratio_side] = find_ratio_opp(...
        ratio_0, ...
        subSubLengthX1{iSub}, subSubLengthY1{iSub}, ...
           minSubSubX1{iSub},    maxSubSubX1{iSub}, ...
           minSubSubY1{iSub},    maxSubSubY1{iSub}, ...
         integrations1{iSub},   surfaceArea1{iSub});    
    ratio_b = [ratio_corner, ratio_0]; 
    ratio_c = [ratio_0, ratio_side];
    
    % Find residuals of ratio_b and ratio_c
    res_b = calc_res(...
        ratio_b(1), ratio_b(2), ...
        subSubLengthX1{iSub}, subSubLengthY1{iSub}, ...
           minSubSubX1{iSub},    maxSubSubX1{iSub}, ...
           minSubSubY1{iSub},    maxSubSubY1{iSub}, ...
         integrations1{iSub},   surfaceArea1{iSub});
    res_c = calc_res(...
        ratio_c(1), ratio_c(2), ...
        subSubLengthX1{iSub}, subSubLengthY1{iSub}, ...
           minSubSubX1{iSub},    maxSubSubX1{iSub}, ...
           minSubSubY1{iSub},    maxSubSubY1{iSub}, ...
         integrations1{iSub},   surfaceArea1{iSub});
    
    % Start triangular bisection
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts
        
        % Calculate next ratio_a and its residual
        ratio_a_next = [...
            (ratio_b(1) + ratio_c(1)) / 2, ...
            (ratio_b(2) + ratio_c(2)) / 2];
        
        res_a_next = calc_res(...
            ratio_a_next(1), ratio_a_next(2), ...
            subSubLengthX1{iSub}, subSubLengthY1{iSub}, ...
               minSubSubX1{iSub},    maxSubSubX1{iSub}, ...
               minSubSubY1{iSub},    maxSubSubY1{iSub}, ...
             integrations1{iSub},   surfaceArea1{iSub});
        
        % Check for convergence
        if abs(ratio_a_next(1) - ratio_a(1)) / ratio_a(1) < areaEps && ...
            abs(ratio_a_next(2) - ratio_a(2)) / ratio_a(2) < areaEps
            
            converged = 1;
            
        end               
        
        % Remove point with same sign as ratio_a_next
        if res_a_next / res_b < 0 && res_a_next / res_c > 0
            
            % ratio_a becomes new ratio_c
            ratio_c = ratio_a;
            
            res_c = res_a;
            
            
        elseif res_a_next / res_c < 0 && res_a_next / res_b > 0
            
            % ratio_a becomes new ratio_b
            ratio_b = ratio_a;
            
            res_b = res_a;
            
            % If signs are the same on opposite side
        else
            
            % Calculate distances of ratio_b and ratio_c from the
            % ratio_corner = ratio_side line            
            dist_b = abs(ratio_b(1) - ratio_b(2)) / sqrt(2);
            
            dist_c = abs(ratio_c(1) - ratio_c(2)) / sqrt(2);
            
            % If distance of ratio_b furhter from the ratio_corner = ratio_side
            % line then remove ratio_b
            if dist_b > dist_c
                
                % ratio_a becomes new ratio_b
                ratio_b = ratio_a;
                
                res_b = res_a;
                                
                % Remove ratio_c
            elseif dist_c > dist_b
                
                % ratio_a becomes new ratio_c
                ratio_c = ratio_a;
                
                res_c = res_a;
                
                % Ratios are equidistant from the ratio_corner = ratio_side
                % line -- arbitrarily remove ratio_c
            else
                
                % ratio_a becomes new ratio_c
                ratio_c = ratio_a;
                
                res_c = res_a;
                
            end
            
        end
        
        % ratio_a_next becomes new ratio_a
        ratio_a = ratio_a_next;
        
        res_a = res_a_next;
        
        counter = counter + 1;
        
    end

    % Solve quadrature directions and weight using converged ratio    
    xPosTemp = [...
        maxSubSubX1{iSub}(1) - subSubLengthX1{iSub}(1) * ratio_a_next(1), ...
       (maxSubSubX1{iSub}(2) +    minSubSubX1{iSub}(2)) / 2, ...
        minSubSubX1{iSub}(3) + subSubLengthX1{iSub}(3) * ratio_a_next(1), ...
        ...
        maxSubSubX1{iSub}(4) - subSubLengthX1{iSub}(4) * ratio_a_next(2), ...
       (maxSubSubX1{iSub}(5) +    minSubSubX1{iSub}(5)) / 2, ...
        minSubSubX1{iSub}(6) + subSubLengthX1{iSub}(6) * ratio_a_next(2), ...
        ...
        maxSubSubX1{iSub}(7) - subSubLengthX1{iSub}(7) * ratio_a_next(1), ...
       (maxSubSubX1{iSub}(8) +    minSubSubX1{iSub}(8)) / 2, ...
        minSubSubX1{iSub}(9) + subSubLengthX1{iSub}(9) * ratio_a_next(1)];
    
    yPosTemp = [...
        maxSubSubY1{iSub}(1) - subSubLengthY1{iSub}(1) * ratio_a_next(1), ...
        maxSubSubY1{iSub}(2) - subSubLengthY1{iSub}(2) * ratio_a_next(2), ...
        maxSubSubY1{iSub}(3) - subSubLengthY1{iSub}(3) * ratio_a_next(1), ...
        ...
        (maxSubSubY1{iSub}(4) + minSubSubY1{iSub}(4)) / 2, ...
        (maxSubSubY1{iSub}(5) + minSubSubY1{iSub}(5)) / 2, ...
        (maxSubSubY1{iSub}(6) + minSubSubY1{iSub}(6)) / 2, ...
        ...
        minSubSubY1{iSub}(7) + subSubLengthY1{iSub}(7) * ratio_a_next(1), ...
        minSubSubY1{iSub}(8) + subSubLengthY1{iSub}(8) * ratio_a_next(2), ...
        minSubSubY1{iSub}(9) + subSubLengthY1{iSub}(9) * ratio_a_next(1)];
    
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations1{iSub}, constants);
    
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
    
    if converged == 0
        
        error('did not converge');
        
    end    
    
    % Store quadrature data
    xPos1{iSub}     = xPosTemp;
    yPos1{iSub}     = yPosTemp;
    gamma1{iSub, 1} = gammaTemp;
    theta1{iSub, 1} = thetaTemp;
    weights1{iSub}  = weightsTemp;
    totWeights = totWeights + sum(weights1{iSub});    
end

%% Generate quadrature for face 2 (x-z plane)

% Initialize storage
xPos2    = cell(numSubSq2, 1);
yPos2    = cell(numSubSq2, 1);
gamma2   = cell(numSubSq2, 1);
theta2   = cell(numSubSq2, 1);
weights2 = cell(numSubSq2, 1);

% Define iteration parameters
areaEps   = 1e-15;
maxCounts = 10000;

% Go through each sub-square
for iSub = 1 : numSubSq2
        
    % Find ratio_a (equal corner and side ratios) minimizing residual
    [ratio_0, res_a] = find_ratio_eq( ...
        subSubLengthX2{iSub}, subSubLengthY2{iSub}, ...
           minSubSubX2{iSub},    maxSubSubX2{iSub}, ...
           minSubSubY2{iSub},    maxSubSubY2{iSub}, ...
         integrations2{iSub},   surfaceArea2{iSub});
    ratio_a = [ratio_0, ratio_0];
    
    % Find ratio_b and ratio_c having opposite residual signs
    [ratio_corner, ratio_side] = find_ratio_opp(...
        ratio_0, ...
        subSubLengthX2{iSub}, subSubLengthY2{iSub}, ...
           minSubSubX2{iSub},    maxSubSubX2{iSub}, ...
           minSubSubY2{iSub},    maxSubSubY2{iSub}, ...
         integrations2{iSub},   surfaceArea2{iSub});    
    ratio_b = [ratio_corner, ratio_0]; 
    ratio_c = [ratio_0, ratio_side];
    
    % Find residuals of ratio_b and ratio_c
    res_b = calc_res(...
        ratio_b(1), ratio_b(2), ...
        subSubLengthX2{iSub}, subSubLengthY2{iSub}, ...
           minSubSubX2{iSub},    maxSubSubX2{iSub}, ...
           minSubSubY2{iSub},    maxSubSubY2{iSub}, ...
         integrations2{iSub},   surfaceArea2{iSub});
    res_c = calc_res(...
        ratio_c(1), ratio_c(2), ...
        subSubLengthX2{iSub}, subSubLengthY2{iSub}, ...
           minSubSubX2{iSub},    maxSubSubX2{iSub}, ...
           minSubSubY2{iSub},    maxSubSubY2{iSub}, ...
         integrations2{iSub},   surfaceArea2{iSub});
    
    % Start triangular bisection
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts
        
        % Calculate next ratio_a and its residual
        ratio_a_next = [...
            (ratio_b(1) + ratio_c(1)) / 2, ...
            (ratio_b(2) + ratio_c(2)) / 2];
        
        res_a_next = calc_res(...
            ratio_a_next(1), ratio_a_next(2), ...
            subSubLengthX2{iSub}, subSubLengthY2{iSub}, ...
               minSubSubX2{iSub},    maxSubSubX2{iSub}, ...
               minSubSubY2{iSub},    maxSubSubY2{iSub}, ...
             integrations2{iSub},   surfaceArea2{iSub});
        
        % Check for convergence
        if abs(ratio_a_next(1) - ratio_a(1)) / ratio_a(1) < areaEps && ...
            abs(ratio_a_next(2) - ratio_a(2)) / ratio_a(2) < areaEps
            
            converged = 1;
            
        end               
        
        % Remove point with same sign as ratio_a_next
        if res_a_next / res_b < 0 && res_a_next / res_c > 0
            
            % ratio_a becomes new ratio_c
            ratio_c = ratio_a;
            
            res_c = res_a;
            
            
        elseif res_a_next / res_c < 0 && res_a_next / res_b > 0
            
            % ratio_a becomes new ratio_b
            ratio_b = ratio_a;
            
            res_b = res_a;
            
            % If signs are the same on opposite side
        else
            
            % Calculate distances of ratio_b and ratio_c from the
            % ratio_corner = ratio_side line            
            dist_b = abs(ratio_b(1) - ratio_b(2)) / sqrt(2);
            
            dist_c = abs(ratio_c(1) - ratio_c(2)) / sqrt(2);
            
            % If distance of ratio_b furhter from the ratio_corner = ratio_side
            % line then remove ratio_b
            if dist_b > dist_c
                
                % ratio_a becomes new ratio_b
                ratio_b = ratio_a;
                
                res_b = res_a;
                                
                % Remove ratio_c
            elseif dist_c > dist_b
                
                % ratio_a becomes new ratio_c
                ratio_c = ratio_a;
                
                res_c = res_a;
                
                % Ratios are equidistant from the ratio_corner = ratio_side
                % line -- arbitrarily remove ratio_c
            else
                
                % ratio_a becomes new ratio_c
                ratio_c = ratio_a;
                
                res_c = res_a;
                
            end
            
        end
        
        % ratio_a_next becomes new ratio_a
        ratio_a = ratio_a_next;
        
        res_a = res_a_next;
        
        counter = counter + 1;
        
    end

    % Solve quadrature directions and weight using converged ratio    
    xPosTemp = [...
        maxSubSubX2{iSub}(1) - subSubLengthX2{iSub}(1) * ratio_a_next(1), ...
       (maxSubSubX2{iSub}(2) +    minSubSubX2{iSub}(2)) / 2, ...
        minSubSubX2{iSub}(3) + subSubLengthX2{iSub}(3) * ratio_a_next(1), ...
        ...
        maxSubSubX2{iSub}(4) - subSubLengthX2{iSub}(4) * ratio_a_next(2), ...
       (maxSubSubX2{iSub}(5) +    minSubSubX2{iSub}(5)) / 2, ...
        minSubSubX2{iSub}(6) + subSubLengthX2{iSub}(6) * ratio_a_next(2), ...
        ...
        maxSubSubX2{iSub}(7) - subSubLengthX2{iSub}(7) * ratio_a_next(1), ...
       (maxSubSubX2{iSub}(8) +    minSubSubX2{iSub}(8)) / 2, ...
        minSubSubX2{iSub}(9) + subSubLengthX2{iSub}(9) * ratio_a_next(1)];
    
    yPosTemp = [...
        maxSubSubY2{iSub}(1) - subSubLengthY2{iSub}(1) * ratio_a_next(1), ...
        maxSubSubY2{iSub}(2) - subSubLengthY2{iSub}(2) * ratio_a_next(2), ...
        maxSubSubY2{iSub}(3) - subSubLengthY2{iSub}(3) * ratio_a_next(1), ...
        ...
        (maxSubSubY2{iSub}(4) + minSubSubY2{iSub}(4)) / 2, ...
        (maxSubSubY2{iSub}(5) + minSubSubY2{iSub}(5)) / 2, ...
        (maxSubSubY2{iSub}(6) + minSubSubY2{iSub}(6)) / 2, ...
        ...
        minSubSubY2{iSub}(7) + subSubLengthY2{iSub}(7) * ratio_a_next(1), ...
        minSubSubY2{iSub}(8) + subSubLengthY2{iSub}(8) * ratio_a_next(2), ...
        minSubSubY2{iSub}(9) + subSubLengthY2{iSub}(9) * ratio_a_next(1)];
    
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations2{iSub}, constants);
    
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
    
    if converged == 0
        
        error('did not converge');
        
    end    
    
    % Store quadrature data
    xPos2{iSub}     = xPosTemp;
    yPos2{iSub}     = yPosTemp;
    gamma2{iSub, 1} = gammaTemp;
    theta2{iSub, 1} = thetaTemp;
    weights2{iSub}  = weightsTemp;
    totWeights = totWeights + sum(weights2{iSub});    
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
xPos3    = cell(numSubSq3, 1);
yPos3    = cell(numSubSq3, 1);
gamma3   = cell(numSubSq3, 1);
theta3   = cell(numSubSq3, 1);
weights3 = cell(numSubSq3, 1);

% Define iteration parameters
areaEps   = 1e-15;
maxCounts = 10000;

% Go through each sub-square
for iSub = 1 : numSubSq3
        
    % Find ratio_a (equal corner and side ratios) minimizing residual
    [ratio_0, res_a] = find_ratio_eq( ...
        subSubLengthX3{iSub}, subSubLengthY3{iSub}, ...
        minSubSubX3{iSub}, maxSubSubX3{iSub}, ...
        minSubSubY3{iSub}, maxSubSubY3{iSub}, ...
        integrations3{iSub}, surfaceArea3{iSub});
    ratio_a = [ratio_0, ratio_0];
    
    % Find ratio_b and ratio_c having opposite residual signs
    [ratio_corner, ratio_side] = find_ratio_opp(...
        ratio_0, ...
        subSubLengthX3{iSub}, subSubLengthY3{iSub}, ...
        minSubSubX3{iSub}, maxSubSubX3{iSub}, ...
        minSubSubY3{iSub}, maxSubSubY3{iSub}, ...
        integrations3{iSub}, surfaceArea3{iSub});    
    ratio_b = [ratio_corner, ratio_0]; 
    ratio_c = [ratio_0, ratio_side];
    
    % Find residuals of ratio_b and ratio_c
    res_b = calc_res(...
        ratio_b(1), ratio_b(2), ...
        subSubLengthX3{iSub}, subSubLengthY3{iSub}, ...
        minSubSubX3{iSub}, maxSubSubX3{iSub}, ...
        minSubSubY3{iSub}, maxSubSubY3{iSub}, ...
        integrations3{iSub}, surfaceArea3{iSub});
    res_c = calc_res(...
        ratio_c(1), ratio_c(2), ...
        subSubLengthX3{iSub}, subSubLengthY3{iSub}, ...
        minSubSubX3{iSub}, maxSubSubX3{iSub}, ...
        minSubSubY3{iSub}, maxSubSubY3{iSub}, ...
        integrations3{iSub}, surfaceArea3{iSub});
    
    % Start triangular bisection
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts
        
        % Calculate next ratio_a and its residual
        ratio_a_next = [...
            (ratio_b(1) + ratio_c(1)) / 2, ...
            (ratio_b(2) + ratio_c(2)) / 2];
        
        res_a_next = calc_res(...
            ratio_a_next(1), ratio_a_next(2), ...
            subSubLengthX3{iSub}, subSubLengthY3{iSub}, ...
            minSubSubX3{iSub}, maxSubSubX3{iSub}, ...
            minSubSubY3{iSub}, maxSubSubY3{iSub}, ...
            integrations3{iSub}, surfaceArea3{iSub});
        
        % Check for convergence
        if abs(ratio_a_next(1) - ratio_a(1)) / ratio_a(1) < areaEps && ...
            abs(ratio_a_next(2) - ratio_a(2)) / ratio_a(2) < areaEps
            
            converged = 1;
            
        end               
        
        % Remove point with same sign as ratio_a_next
        if res_a_next / res_b < 0 && res_a_next / res_c > 0
            
            % ratio_a becomes new ratio_c
            ratio_c = ratio_a;
            
            res_c = res_a;
            
            
        elseif res_a_next / res_c < 0 && res_a_next / res_b > 0
            
            % ratio_a becomes new ratio_b
            ratio_b = ratio_a;
            
            res_b = res_a;
            
            % If signs are the same on opposite side
        else
            
            % Calculate distances of ratio_b and ratio_c from the
            % ratio_corner = ratio_side line            
            dist_b = abs(ratio_b(1) - ratio_b(2)) / sqrt(2);
            
            dist_c = abs(ratio_c(1) - ratio_c(2)) / sqrt(2);
            
            % If distance of ratio_b furhter from the ratio_corner = ratio_side
            % line then remove ratio_b
            if dist_b > dist_c
                
                % ratio_a becomes new ratio_b
                ratio_b = ratio_a;
                
                res_b = res_a;
                                
                % Remove ratio_c
            elseif dist_c > dist_b
                
                % ratio_a becomes new ratio_c
                ratio_c = ratio_a;
                
                res_c = res_a;
                
                % Ratios are equidistant from the ratio_corner = ratio_side
                % line -- arbitrarily remove ratio_c
            else
                
                % ratio_a becomes new ratio_c
                ratio_c = ratio_a;
                
                res_c = res_a;
                
            end
            
        end
        
        % ratio_a_next becomes new ratio_a
        ratio_a = ratio_a_next;
        
        res_a = res_a_next;
        
        counter = counter + 1;
        
    end

    % Solve quadrature directions and weight using converged ratio    
    xPosTemp = [...
        maxSubSubX3{iSub}(1) - subSubLengthX3{iSub}(1) * ratio_a_next(1), ...
        (maxSubSubX3{iSub}(2) + minSubSubX3{iSub}(2)) / 2, ...
        minSubSubX3{iSub}(3) + subSubLengthX3{iSub}(3) * ratio_a_next(1), ...
        ...
        maxSubSubX3{iSub}(4) - subSubLengthX3{iSub}(4) * ratio_a_next(2), ...
        (maxSubSubX3{iSub}(5) + minSubSubX3{iSub}(5)) / 2, ...
        minSubSubX3{iSub}(6) + subSubLengthX3{iSub}(6) * ratio_a_next(2), ...
        ...
        maxSubSubX3{iSub}(7) - subSubLengthX3{iSub}(7) * ratio_a_next(1), ...
        (maxSubSubX3{iSub}(8) + minSubSubX3{iSub}(8)) / 2, ...
        minSubSubX3{iSub}(9) + subSubLengthX3{iSub}(9) * ratio_a_next(1)];
    
    yPosTemp = [...
        maxSubSubY3{iSub}(1) - subSubLengthY3{iSub}(1) * ratio_a_next(1), ...
        maxSubSubY3{iSub}(2) - subSubLengthY3{iSub}(2) * ratio_a_next(2), ...
        maxSubSubY3{iSub}(3) - subSubLengthY3{iSub}(3) * ratio_a_next(1), ...
        ...
        (maxSubSubY3{iSub}(4) + minSubSubY3{iSub}(4)) / 2, ...
        (maxSubSubY3{iSub}(5) + minSubSubY3{iSub}(5)) / 2, ...
        (maxSubSubY3{iSub}(6) + minSubSubY3{iSub}(6)) / 2, ...
        ...
        minSubSubY3{iSub}(7) + subSubLengthY3{iSub}(7) * ratio_a_next(1), ...
        minSubSubY3{iSub}(8) + subSubLengthY3{iSub}(8) * ratio_a_next(2), ...
        minSubSubY3{iSub}(9) + subSubLengthY3{iSub}(9) * ratio_a_next(1)];
    
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations3{iSub}, constants);
    
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
    
    if converged == 0
        
        error('did not converge');
        
    end    
    
    % Store quadrature data
    xPos3{iSub}     = xPosTemp;
    yPos3{iSub}     = yPosTemp;
    gamma3{iSub, 1} = gammaTemp;
    theta3{iSub, 1} = thetaTemp;
    weights3{iSub}  = weightsTemp;
    totWeights = totWeights + sum(weights3{iSub});    
end

% Rotate to face 3
for i = 1 : numSubSq3
    for j = 1 : 9
        gamma3{i}(j) = atan(yPos3{i}(j) / xPos3{i}(j));
        theta3{i}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos3{i}(j) ^ 2 + xPos3{i}(j) ^ 2)));
    end
end

%% Maximum surface area error
test = zeros(numSubSq1 + numSubSq2 + numSubSq3, 1);
counter = 1;
for i = 1 : numSubSq1
    test(counter) = abs(weights1{i}(5) - surfaceArea1{i}(5)) / surfaceArea1{i}(5);
    counter = counter + 1;
end
for i = 1 : numSubSq2
    test(counter) = abs(weights2{i}(5) - surfaceArea2{i}(5)) / surfaceArea2{i}(5);
    counter = counter + 1;
end
for i = 1 : numSubSq3
    test(counter) = abs(weights3{i}(5) - surfaceArea3{i}(5)) / surfaceArea3{i}(5);
    counter = counter + 1;
end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

%% Store QDFE-ratio quadrature data
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