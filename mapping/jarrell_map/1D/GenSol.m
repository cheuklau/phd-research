function sol = GenSol(quad, functionNumber)

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
    
    % Smooth linear function
elseif functionNumber == 3
    for i = 1 : numRegions
        for j = 1 : 2
            sol{i}(j) = 1 + mu{i}(j);
        end
    end
    
    % Exponential function with discontinuity at mu=0
elseif functionNumber == 4
    for i = 1 : numRegions
        for j = 1 : 2
            sol{i}(j) = exp(-1 / mu{i}(j));
        end
    end
    
end


end