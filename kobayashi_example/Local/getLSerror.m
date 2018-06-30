function [num_error, num_sol, meshSize] = getLSerror(ref_LS, refSol, prob)

% Add path
path(path, './Pregen Functions');

fprintf('Starting LS error calculation \n');

% Number of directions for each refinement level
for i = 1 : size(ref_LS, 2)
   
    numLines(i) = ref_LS(i) * (ref_LS(i) + 2) / 8;
    
end

% Initialize relative error
intError = [];

% Import LS data
dataLS = dlmread('LS.DAT');

% Go through each refinement level
for i = 1 : size(ref_LS, 2);
    
    fprintf('Refinement level: %i \n', ref_LS(i));
    
    % Starting line
    lineStart = 1;
    for j = 1 : i - 1
        lineStart = lineStart + numLines(j);
    end
    
    % Store data values    
    omegaX      = zeros(1, numLines(i));
    omegaY      = zeros(1, numLines(i));
    omegaZ      = zeros(1, numLines(i));
    weights     = zeros(1, numLines(i));
    weightTotal = 0;
    for j = 1 : numLines(i)
        lineUse = lineStart + j - 1;
        omegaX(j)   = dataLS(lineUse, 1);
        omegaY(j)   = dataLS(lineUse, 2);
        omegaZ(j)   = dataLS(lineUse, 3);
        weights(j)  = dataLS(lineUse, 4);
        weightTotal = weightTotal + weights(j);
    end
    
    % Normalize weights
    for j = 1 : numLines(i)
       weights(j) = weights(j) * (4 * pi / 8) / weightTotal; 
    end
    
    % Store quadrature
    quad.omegaX  = omegaX;
    quad.omegaY  = omegaY;
    quad.omegaZ  = omegaZ;
    quad.weights = weights;
    
    % Calculate error
    [num_error(i) num_sol(i)] = SphInt(quad, refSol, prob);
    
    % Angular mesh size
    meshSize(i) = 1 / sqrt(8 * numLines(i));
        
end

rmpath('./Pregen Functions');

end