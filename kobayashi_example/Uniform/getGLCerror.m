function [num_error, num_sol, meshSize] = getGLCerror(ref_GLC, refSol, prob)

% Add path
path(path, './Pregen Functions');

fprintf('Starting GLC error calculation \n');

% Number of directions for each refinement level
for i = 1 : size(ref_GLC, 2)
   
    numLines(i) = ref_GLC(i) * (ref_GLC(i) + 1) / 2;
    
end

% Initialize relative error
intError = [];

% Import GLC data
dataGLC = dlmread('GLC.DAT');

% Go through each refinement level
for i = 1 : size(ref_GLC, 2);
    
    fprintf('Refinement level: %i \n', ref_GLC(i));
    
    % Starting line
    lineStart = 1;
    for j = 1 : i - 1
        lineStart = lineStart + numLines(j);
    end
    
    % Store data values
    omegaX = zeros(1, numLines(i));
    omegaY = zeros(1, numLines(i));
    omegaZ = zeros(1, numLines(i));
    weights = zeros(1, numLines(i));
    weightTotal = 0;
    for j = 1 : numLines(i)
        lineUse = lineStart + j - 1;
        omegaX(j) = dataGLC(lineUse, 1);
        omegaY(j) = dataGLC(lineUse, 2);
        omegaZ(j) = dataGLC(lineUse, 3);
        weights(j) = dataGLC(lineUse, 4);
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