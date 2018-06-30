% Generate QDFE-ratio quadrature
function [rho_corner, rho_side] = InitialRatio(squareInfo)

% Retrieve needed square properties
minSubSubX   = squareInfo.minSubSubX;
maxSubSubX   = squareInfo.maxSubSubX;
minSubSubY   = squareInfo.minSubSubY;
maxSubSubY   = squareInfo.maxSubSubY;
subSubLengthX = squareInfo.subSubLengthX;
subSubLengthY = squareInfo.subSubLengthY;
surfaceArea  = squareInfo.surfaceArea;
integrations = squareInfo.integrations;
numSubSq     = squareInfo.numSubSq;

% Initialize storage
rho_corner = zeros(numSubSq, 1);
rho_side   = zeros(numSubSq, 1);

% Define iteration parameters
areaEps   = 1e-15;
maxCounts = 1000;

% Go through each sub-square
for iSub = 1 : numSubSq
    
    fprintf('Currenlty on sub-square %i of %i \n', iSub, numSubSq);
        
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

    % Store converged ratio
    rho_corner(iSub) = ratio_a_next(1);
    rho_side(iSub)   = ratio_a_next(2);
        
    % Check if ratio is valid
    if min(ratio_a_next) < 0 || max(ratio_a_next) >= 1
        
        error('bad ratio');
    
    end
       
end

end