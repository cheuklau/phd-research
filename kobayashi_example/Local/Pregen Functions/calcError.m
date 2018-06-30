function [num_err num_sol] = calcError(quad, refSol, prob)

%% Problem information
x = prob.x;
y = prob.y;
z = prob.z;
P0 = [x y z];
sigma_a_source = prob.sigma_a_source;
sigma_a_box    = prob.sigma_a_box;
sigma_a_duct   = prob.sigma_a_duct;
source         = prob.source;

%% Gather quadrature information
omegaX  = quad.omegaX;
omegaY  = quad.omegaY;
omegaZ  = quad.omegaZ;
weights = quad.weights;
numDirs = size(omegaX, 2);

%% Numerical solution using quadrature
num_sol = 0;

% Go through each direction
for i = 1 : numDirs
    
    %% Initialize variables
    
    % Find second point on line
    x_long = x - omegaX(i);
    y_long = y - omegaY(i);
    z_long = z - omegaZ(i);
    P1 = [x_long, y_long, z_long];
    
    % Initialize pathlengths
    p_reg2   = 0;
    p_reg3   = 0;
    p_reg4   = 0;
    p_source = 0;
    
    %% Find pathlength through region #1
    
    % Set bounds
    x_bounds = [40 50];
    y_bounds = [60 110];
    z_bounds = [40 50];
    
    % Set min and max vertices
    Vmin = [x_bounds(1) y_bounds(1) z_bounds(1)];
    
    % Calculate pathlength to exit
    p_reg1 = calc_exit(Vmin, P0, P1, x_bounds, y_bounds, z_bounds, P0, i);
    
    
    %% Find pathlength through region #2
    
    % Set bounds
    x_bounds = [40 50];
    y_bounds = [60 70];
    z_bounds = [0 40];
    
    % Set min and max vertices
    Vmin = [x_bounds(1) y_bounds(1) z_bounds(1)];
    Vmax = [x_bounds(2) y_bounds(2) z_bounds(2)];
    
    % Intersection point for entering face
    p_enter = calc_enter(Vmax, P0, P1, x_bounds, y_bounds, z_bounds);
    
    if max(p_enter) > 0 
        
        p_reg2 = calc_exit(Vmin, P0, P1, x_bounds, y_bounds, z_bounds, p_enter, i);
        
    end
    
    %% Find pathlength through region #3
    
    % Set bounds
    x_bounds = [0  40];
    y_bounds = [60 70];
    z_bounds = [0  20];
    
    % Set min and max vertices
    Vmin = [x_bounds(1) y_bounds(1) z_bounds(1)];
    Vmax = [x_bounds(2) y_bounds(2) z_bounds(2)];
    
    % Intersection point for entering face
    p_enter = calc_enter(Vmax, P0, P1, x_bounds, y_bounds, z_bounds);
    
    if max(p_enter) > 0 
        
        p_reg3 = calc_exit(Vmin, P0, P1, x_bounds, y_bounds, z_bounds, p_enter, i);
        
    end
    
    % Set bounds
    x_bounds = [0  20];
    y_bounds = [20 60];
    z_bounds = [0  20];
    
    % Set min and max vertices
    Vmin = [x_bounds(1) y_bounds(1) z_bounds(1)];
    Vmax = [x_bounds(2) y_bounds(2) z_bounds(2)];
    
    % Intersection point for entering face
    p_enter = calc_enter(Vmax, P0, P1, x_bounds, y_bounds, z_bounds);
    
    if max(p_enter) > 0 
        
        p_reg4 = calc_exit(Vmin, P0, P1, x_bounds, y_bounds, z_bounds, p_enter, i);
        
    end
    
    %% Find pathlength through source
    
    % Set bounds
    x_bounds = [0 20];
    y_bounds = [0 20];
    z_bounds = [0 20];
    
    % Set min and max vertices
    Vmin = [x_bounds(1) y_bounds(1) z_bounds(1)];
    Vmax = [x_bounds(2) y_bounds(2) z_bounds(2)];
    
    % Intersection point for entering face
    p_enter = calc_enter(Vmax, P0, P1, x_bounds, y_bounds, z_bounds);
    
    if max(p_enter) > 0 
        
        p_source = calc_exit(Vmin, P0, P1, x_bounds, y_bounds, z_bounds, p_enter, i);
        
    end
    
    %% Calculate contribution to scalar flux
    if p_source > 0
                
        p_duct = p_reg1 + p_reg2 + p_reg3 + p_reg4;
        
        p_box = sqrt((x - p_enter(1)) ^2 + (y - p_enter(2)) ^ 2 + (z - p_enter(3)) ^ 2) - p_duct;
        
        num_sol = num_sol + ...
            (source / (4 * pi * sigma_a_source)) * (1 - exp(-1 * sigma_a_source * p_source)) * ...
            exp(-1 * sigma_a_box * p_box) * exp(-1 * sigma_a_duct * p_duct) * weights(i);
        
    end
    
    
end

num_err = abs(refSol - num_sol) / refSol;

end