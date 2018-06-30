% Generate QDFE-ratio quadrature
function quadrature = QDFEratioTBS(squareInfo)

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
xPos    = cell(numSubSq, 1);
yPos    = cell(numSubSq, 1);
gamma   = cell(numSubSq, 3);
theta   = cell(numSubSq, 3);
weights = cell(numSubSq, 1);

% Define iteration parameters
areaEps   = 1e-15;
maxCounts = 1000;

% Go through each sub-square
for iSub = 1 : numSubSq
    
    fprintf('Currently on sub-square %i of %i \n', iSub, numSubSq);
    
    % Find ratio_a (equal corner and side ratios) minimizing residual
    [ratio_0, res_a] = find_ratio_eq( ...
        subSubLengthX{iSub}, subSubLengthY{iSub}, ...
        minSubSubX{iSub}, maxSubSubX{iSub}, ...
        minSubSubY{iSub}, maxSubSubY{iSub}, ...
        integrations{iSub}, surfaceArea{iSub});
    ratio_a = [ratio_0, ratio_0];
    
    % Find ratio_b and ratio_c having opposite residual signs
    [ratio_corner, ratio_side] = find_ratio_opp(...
        ratio_0, ...
        subSubLengthX{iSub}, subSubLengthY{iSub}, ...
        minSubSubX{iSub}, maxSubSubX{iSub}, ...
        minSubSubY{iSub}, maxSubSubY{iSub}, ...
        integrations{iSub}, surfaceArea{iSub});    
    ratio_b = [ratio_corner, ratio_0]; 
    ratio_c = [ratio_0, ratio_side];
    
    % Find residuals of ratio_b and ratio_c
    res_b = calc_res(...
        ratio_b(1), ratio_b(2), ...
        subSubLengthX{iSub}, subSubLengthY{iSub}, ...
        minSubSubX{iSub}, maxSubSubX{iSub}, ...
        minSubSubY{iSub}, maxSubSubY{iSub}, ...
        integrations{iSub}, surfaceArea{iSub});
    res_c = calc_res(...
        ratio_c(1), ratio_c(2), ...
        subSubLengthX{iSub}, subSubLengthY{iSub}, ...
        minSubSubX{iSub}, maxSubSubX{iSub}, ...
        minSubSubY{iSub}, maxSubSubY{iSub}, ...
        integrations{iSub}, surfaceArea{iSub});
    
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
            subSubLengthX{iSub}, subSubLengthY{iSub}, ...
            minSubSubX{iSub}, maxSubSubX{iSub}, ...
            minSubSubY{iSub}, maxSubSubY{iSub}, ...
            integrations{iSub}, surfaceArea{iSub});
        
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
        maxSubSubX{iSub}(1) - subSubLengthX{iSub}(1) * ratio_a_next(1), ...
        (maxSubSubX{iSub}(2) + minSubSubX{iSub}(2)) / 2, ...
        minSubSubX{iSub}(3) + subSubLengthX{iSub}(3) * ratio_a_next(1), ...
        ...
        maxSubSubX{iSub}(4) - subSubLengthX{iSub}(4) * ratio_a_next(2), ...
        (maxSubSubX{iSub}(5) + minSubSubX{iSub}(5)) / 2, ...
        minSubSubX{iSub}(6) + subSubLengthX{iSub}(6) * ratio_a_next(2), ...
        ...
        maxSubSubX{iSub}(7) - subSubLengthX{iSub}(7) * ratio_a_next(1), ...
        (maxSubSubX{iSub}(8) + minSubSubX{iSub}(8)) / 2, ...
        minSubSubX{iSub}(9) + subSubLengthX{iSub}(9) * ratio_a_next(1)];
    
    yPosTemp = [...
        maxSubSubY{iSub}(1) - subSubLengthY{iSub}(1) * ratio_a_next(1), ...
        maxSubSubY{iSub}(2) - subSubLengthY{iSub}(2) * ratio_a_next(2), ...
        maxSubSubY{iSub}(3) - subSubLengthY{iSub}(3) * ratio_a_next(1), ...
        ...
        (maxSubSubY{iSub}(4) + minSubSubY{iSub}(4)) / 2, ...
        (maxSubSubY{iSub}(5) + minSubSubY{iSub}(5)) / 2, ...
        (maxSubSubY{iSub}(6) + minSubSubY{iSub}(6)) / 2, ...
        ...
        minSubSubY{iSub}(7) + subSubLengthY{iSub}(7) * ratio_a_next(1), ...
        minSubSubY{iSub}(8) + subSubLengthY{iSub}(8) * ratio_a_next(2), ...
        minSubSubY{iSub}(9) + subSubLengthY{iSub}(9) * ratio_a_next(1)];
    
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations{iSub}, constants);
    
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
    
    if converged == 0
        
        error('did not converge');
        
    end    
    
    % Store quadrature data
    xPos{iSub}     = xPosTemp;
    yPos{iSub}     = yPosTemp;
    gamma{iSub, 1} = gammaTemp;
    theta{iSub, 1} = thetaTemp;
    weights{iSub}  = weightsTemp;
    
end

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
    test(counter) = abs(weights{i}(5) - surfaceArea{i}(5)) / surfaceArea{i}(5);
    counter = counter + 1;
end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

% Rotate to other three faces
for i = 1 : numSubSq
    for j = 1 : 9
        gamma{i, 2}(j) = (pi / 2) - atan(sqrt(3) * yPos{i}(j));
        theta{i, 2}(j) = (pi / 2) - atan(sqrt(3) * xPos{i}(j) / sqrt(1 + 3 * yPos{i}(j) ^ 2));
        gamma{i, 3}(j) = atan(yPos{i}(j) / xPos{i}(j));
        theta{i, 3}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos{i}(j) ^ 2 + xPos{i}(j) ^ 2)));
    end
end

% Store QDFE-ratio quadrature data
quadrature.gamma   = gamma;
quadrature.theta   = theta;
quadrature.weights = weights;

end