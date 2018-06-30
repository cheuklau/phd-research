function [sol, limits] = GenSol(quad, functionNumber)

% Quadrature directions
mu         = quad.mu;
numRegions = size(mu, 2);

% Initialize storage
sol = cell(numRegions, 1);

% Step function with discontinuity within an LDFE region
if functionNumber == 1
    for i = 1 : numRegions
        for j = 1 : 2
            if mu{i}(j) <= -1 / 3
                sol{i}(j) = 1;
            else
                sol{i}(j) = 0;
            end
        end        
    end
    maxTmp = cellfun(@(x) max(x(:)), sol);
    minTmp = cellfun(@(x) min(x(:)), sol);
    psiMax = max(maxTmp);
    psiMin = min(minTmp);

    % Step function with discontinuity at mu=0
elseif functionNumber == 2
    for i = 1 : numRegions
        for j = 1 : 2
            if mu{i}(j) <= 0
                sol{i}(j) = 1;
            else
                sol{i}(j) = 0;
            end
        end   
    end
    maxTmp = cellfun(@(x) max(x(:)), sol);
    minTmp = cellfun(@(x) min(x(:)), sol);
    psiMax = max(maxTmp);
    psiMin = min(minTmp);
    
    % Smooth linear function
elseif functionNumber == 3
    for i = 1 : numRegions
        for j = 1 : 2
            sol{i}(j) = 1 + mu{i}(j);
        end 
    end
    maxTmp = cellfun(@(x) max(x(:)), sol);
    minTmp = cellfun(@(x) min(x(:)), sol);
    psiMax = max(maxTmp);
    psiMin = min(minTmp);
    
    % Exponential function with discontinuity at mu=0
elseif functionNumber == 4
    for i = 1 : numRegions
        for j = 1 : 2
            sol{i}(j) = exp(-1 / mu{i}(j));
        end     
    end
    maxTmp = cellfun(@(x) max(x(:)), sol);
    minTmp = cellfun(@(x) min(x(:)), sol);
    psiMax = max(maxTmp);
    psiMin = min(minTmp);
    
    % Preliminary exam problem
elseif functionNumber == 5
    
    % Problem parameters
    S       = 1.0; % Source strength
    sigma_t = 0.1; % Total cross section
    delta_x = 0.001; % Source thickness
    
    % Analytic solution
    for i = 1 : numRegions
        for j = 1 : 2
           sol{i}(j) = (S / sigma_t) * ...
               (1 - exp(-1 * sigma_t * delta_x / (2 * abs(mu{i}(j)))));
        end         
    end 
    psiMax = S / sigma_t;
    psiMin = 0;
    
end

% Store limits
limits.psiMax = psiMax;
limits.psiMin = psiMin;

end