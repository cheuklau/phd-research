%% Prime-Dual Infeasible Interior-Point Algorithm
% Taken from Numerical Recipes (version 3)
% Input:
% a = coefficient matrix for the constrants
% b = right-hand side of the constraints
% c = coefficients of the objective function to be minimized
% Note: c should be padded with the slack variables

function [status, iter, x] = intpt(a, b, c)

%% Define parameters
MAXITS = 200;        % Maximum number of iterations
TOL    = 1e-6;      % Tolerance for optimality and feasibility
SIGMA  = 0.9;      % Stepsize reduction factor (conservative choice)
DELTA  = 0.02;  % Factor to set centrality parameter mu
BIG    = realmax;    % Largest Matlab value
m      = size(a, 1); % Number of constraints
n      = size(a, 2); % Number of degrees of freedom
status = -1;         % Initialize return status
d      = zeros(n);   % Initialize diagonal matrix X*Z^-1

%% Transpose of a-matrix
at = transpose(a);

%% Factors for convergence test
rpfact = 1 + sqrt(dot(b, b));
rdfact = 1 + sqrt(dot(c, c));

%% Initial point
% Note: Can be outside of feasible region (as long as non-negative)
x = 10000 * ones(n, 1);
z = 10000 * ones(n, 1);
y = 10000 * ones(m, 1);

%% Initialize residual vector norms
normrp_old = BIG;
normrd_old = BIG;

%% Main loop
for iter = 1 : MAXITS
    
    % Primal residual and norm
    rp     = a * x - b;
    normrp = sqrt(dot(rp, rp)) / rpfact;
    
    % Dual residual and norm
    rd     = at * y + z - c;
    normrd = sqrt(dot(rd, rd)) / rdfact;
    
    % Duality gap
    gamma = dot(x, z);
    
    % Duality measure
    mu = DELTA * gamma / n;
    
    % Primal objective
    primal_obj = dot(c, x);
    
    % Duality gap convergence criteria
    gamma_norm = gamma / (1 + abs(primal_obj));
    
    % Optimal solution found
    if normrp < TOL && normrd < TOL && gamma_norm < TOL      
        
        status = 0;
        
        break;
        
        % Primal infeasible
    elseif normrp > 1000 * normrp_old && normrp > TOL        
        
        status = 1;
        
        break;
        
        % Dual infeasible
    elseif normrd > 1000 * normrd_old && normrd > TOL        
        
        status = 2;
        
        break;        
        
    end
    
    % Diagonal matrix entries
    for i = 1 : n          
        
        d(i, i) = x(i) / z(i);               
        
    end
    
    % Calculate A * (X / Z) * At
    adat = a * d * at;
    
    % Form right-hand side
    tempn = x - mu ./ z - d * rd;
    tempm = a * tempn;
    rhs = -1 * rp + tempm;
    
    % Solve for dy
    dy = adat \ rhs;
    
    % Solve for dz
    tempn = at * dy;    
    dz = -1 * tempn - rd;
    
    % Solve for dx
    dx = -1 * d * dz + mu ./ z - x;
    
    % Find factor for primal step size
    alpha_p = 1;
    
    for j = 1 : n        
        
        if x(j) + alpha_p * dx(j) < 0            
            
            alpha_p = -1 * x(j) / dx(j);            
            
        end 
        
    end 
    
    % Reduce alpha by safety factor to ensure non-negativity
    alpha_p = min(alpha_p * SIGMA, 1);
    
    % Find factor for dual step size
    alpha_d = 1;    
    
    for j = 1 : n        
        
        if z(j) + alpha_d * dz(j) < 0            
            
            alpha_d = -1 * z(j) / dz(j);            
            
        end 
        
    end
    
    % Reduce alpha by safety factor to ensure non-negativity    
    alpha_d = min(alpha_d * SIGMA, 1);
    
    % Step to new point
    x = x + alpha_p * dx;
    z = z + alpha_d * dz;
    y = y + alpha_d * dy;
    
    % Update residual norm
    normrp_old = normrp;    
    normrd_old = normrd;
    
end

%% Check for non-convergence
if status == -1
    
    status = 3;
    
end

end