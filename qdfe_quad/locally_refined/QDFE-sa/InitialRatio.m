% Generate QDFE-ratio quadrature
function [rho_corner1, rho_side1, ...
          rho_corner2, rho_side2, ...
          rho_corner3, rho_side3] = InitialRatio(squareInfo)

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
rho_corner1 = zeros(numSubSq1, 1);
rho_side1   = zeros(numSubSq1, 1);

% Define iteration parameters
areaEps   = 1e-15;
maxCounts = 1000;

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

    % Store converged ratio
    rho_corner1(iSub) = ratio_a_next(1);
    rho_side1(iSub) = ratio_a_next(2);
    
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
        
end

%% Generate quadrature for face 2 (x-z plane)

% Initialize storage
rho_corner2 = zeros(numSubSq2, 1);
rho_side2   = zeros(numSubSq2, 1);

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

    % Store converged ratio
    rho_corner2(iSub) = ratio_a_next(1);
    rho_side2(iSub)   = ratio_a_next(2);
    
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
    
end

%% Generate quadrature for face 3 (y-z plane)

% Initialize storage
rho_corner3 = zeros(numSubSq3, 1);
rho_side3   = zeros(numSubSq3, 1);

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

    % Store converged ratio
    rho_corner3(iSub) = ratio_a_next(1);
    rho_side3(iSub)   = ratio_a_next(2);
    
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
    
end

end