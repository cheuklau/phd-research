% Find side and corner ratios with opposite signs for initial guess
function [ratio_corner, ratio_side] = find_ratio_opp(...
    ratio_eq, ...
    subSubLengthX, subSubLengthY, ...
    minSubSubX, maxSubSubX, ...
    minSubSubY, maxSubSubY, ...
    integrations, surfaceArea)

% Amount to increment ratios
delta = 1e-2;

% Start search for side and corner ratios
found = 0;
while found == 0 && ratio_eq + delta <= 1.0
    
    % New corner ratio and residual
    ratio_corner = ratio_eq + delta;
    res_corner = calc_res(...
        ratio_corner, ratio_eq, ...
        subSubLengthX, subSubLengthY, ...
        minSubSubX, maxSubSubX, ...
        minSubSubY, maxSubSubY, ...
        integrations, surfaceArea);
    
    % New side ratio and residual 
    ratio_side = ratio_eq - delta;
    res_side = calc_res(...
        ratio_eq, ratio_side, ...
        subSubLengthX, subSubLengthY, ...
        minSubSubX, maxSubSubX, ...
        minSubSubY, maxSubSubY, ...
        integrations, surfaceArea);
    
    % Determine if the residal signs are different
    if res_corner / res_side < 0
        
        found = 1;
        
        % Otherwise increase change in ratio
    else
        
        delta = delta + 1e-2;
        
    end
        
end

% Check if opposite signs could not be found
if found == 0
    
    error('could not find opposite sign ratios!');
    
end

end