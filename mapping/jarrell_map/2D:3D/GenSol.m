function sol = GenSol(quad, functionNumber, problem)

% Number of directions
numDirs = size(quad, 1);

% Initialize storage
sol = zeros(numDirs, 1);

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
        
        quad_mag = sqrt(quad{i}(1) ^ 2 + quad{i}(2) ^ 2 + quad{i}(3) ^ 2);
        
        t_i = acos(dot_prod / (d * quad_mag));
        
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
    
end

end