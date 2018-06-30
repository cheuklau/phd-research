function [sol, limits] = GenSol(quad, functionNumber, problem)

%% Number of directions
numDirs = size(quad, 1);

%% Initialize storage
sol = zeros(numDirs, 1);

%% Calculate reference solution for given quadrature

% Linear function
if functionNumber == 1
    
    for i = 1 : numDirs
        
        sol(i) = 1 + quad{i}(1) + quad{i}(2) + quad{i}(3);
        
    end

    % Quadratic function
elseif functionNumber == 2
    
    for i = 1 : numDirs
        
        sol(i) = 1 + quad{i}(1) ^ 5 + quad{i}(2) ^ 5 + quad{i}(3) ^ 5;
        
    end
    
    % Spherical source function
elseif functionNumber == 3
    
    % Retrieve problem information
    R = problem.source_radius;
    Q = problem.source_strength;
    S = problem.source_sigma_t;
    P = problem.detector_pos;
    
    % Distance from detector to origin
    d = sqrt(P(1) ^ 2 + P(2) ^ 2 + P(3) ^ 2);
    
    % Maximum angle between detector and source
    t_max = atan(R / d);
        
    % Go through each quadrature direction
    for i = 1 : numDirs
       
        % Angle between quadrature direction and source
        dot_prod = P(1) * quad{i}(1) + P(2) * quad{i}(2) + P(3) * quad{i}(3);
                
        t_i = acos(dot_prod / d);
        
        % Non-zero angular flux if quadrature direction passes through sphere
        if t_i < t_max
                
            % Calculate chord length
            alpha = pi / 2 - t_i;
            
            r = d / tan(alpha);
            
            beta = asin(sin(alpha) * (r / R));
            
            delta = pi - beta - alpha;
            
            chord_1 = sin(delta) * r / sin(beta);
            
            gamma = pi  - alpha;
            
            xi = asin(r * sin(gamma) / R);
            
            chi = pi - gamma - xi;
            
            chord_2 = sin(chi) * r / sin(xi);
            
            chord = chord_1 + chord_2;
            
            % Calculate angular flux solution
            sol(i) = (Q / S) * (1 - exp(-1 * S * chord));
            
            % Zero angular flux if quadrature flux does not pass through sphere
        else
            
            sol(i) = 0;
            
        end

    end
  
    % Discontinuous function
elseif functionNumber == 4
    
    % Gamma and theta limits
    min_gamma = problem.min_gamma;
    max_gamma = problem.max_gamma;
    min_theta = problem.min_theta;
    max_theta = problem.max_theta;
    
    for i = 1 : numDirs
       
        % Calculate theta and gamma of current direction
        theta = acos(quad{i}(3));
        gamma = acos(quad{i}(1) / sin(theta));
        
        % Determine if current direction is within patch of interest
        if theta > min_theta && theta < max_theta && ...
                gamma > min_gamma && gamma < max_gamma
           
            sol(i) = 1;
            
        else
            
            sol(i) = 0;
            
        end
        
    end
    
end

%% Calculate maximum and minimum solution using a fine mesh
gridSize = 1000;
gamma = linspace(0, pi / 2, gridSize);
theta = linspace(0, pi / 2, gridSize);
omegaXref = zeros(gridSize);
omegaYref = zeros(gridSize);
omegaZref = zeros(gridSize);
solRef = zeros(gridSize);
for i = 1 : gridSize
    for j = 1 : gridSize
        
        % Calculate directional cosines
        omegaXref(i, j) = cos(gamma(i)) * sin(theta(j));
        omegaYref(i, j) = sin(gamma(i)) * sin(theta(j));
        omegaZref(i, j) = cos(theta(j));
        
        % Linear function
        if functionNumber == 1
            
            solRef(i, j) = 1 + omegaXref(i, j) + omegaYref(i, j) + omegaZref(i, j);
            
            % Quadratic function
        elseif functionNumber == 2
            
            solRef(i, j) = 1 + omegaXref(i, j) ^ 5 + omegaYref(i, j) ^ 5 + omegaZref(i, j) ^ 5;
            
            % Spherical source function
        elseif functionNumber == 3
            
            % Retrieve problem information
            R = problem.source_radius;
            Q = problem.source_strength;
            S = problem.source_sigma_t;
            P = problem.detector_pos;
            
            % Distance from detector to origin
            d = sqrt(P(1) ^ 2 + P(2) ^ 2 + P(3) ^ 2);
            
            % Maximum angle between detector and source
            t_max = atan(R / d);
            
            % Angle between quadrature direction and source
            dot_prod = ...
                P(1) * omegaXref(i, j) + ...
                P(2) * omegaYref(i, j) + ...
                P(3) * omegaZref(i, j);          
            
            t_i = acos( dot_prod / d);
            
            % Non-zero angular flux if quadrature direction passes through sphere
            if t_i < t_max
                
                % Calculate chord length
                alpha = pi / 2 - t_i;
                
                r = d / tan(alpha);
                
                beta = asin(sin(alpha) * (r / R));
                
                delta = pi - beta - alpha;
                
                chord_1 = sin(delta) * r / sin(beta);
                
                gamma_2 = pi  - alpha;
                
                xi = asin(r * sin(gamma_2) / R);
                
                chi = pi - gamma_2 - xi;
                
                chord_2 = sin(chi) * r / sin(xi);
                
                chord = chord_1 + chord_2;
                
                % Calculate angular flux solution
                solRef(i, j) = (Q / S) * (1 - exp(-1 * S * chord));
                
                % Zero angular flux if quadrature flux does not pass through sphere
            else
                
                solRef(i, j) = 0;
                
            end
            
        end                
        
    end
    
end

% Store limits
if functionNumber ~= 4
    
    limits.psiMax = max(max(solRef));
    limits.psiMin = min(min(solRef));

else
    
    limits.psiMax = 1;
    limits.psiMin = 0;
    
end


end