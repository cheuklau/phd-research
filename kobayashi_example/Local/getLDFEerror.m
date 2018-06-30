function [num_error, num_sol, meshSize] = getLDFEerror(ref_LDFE, refSol, prob)

%% Add path
path(path, './Pregen Functions');

fprintf('Starting LDFE error calculation \n');

%% Number of directions for each refinement level
numLines = ref_LDFE;

%% Initialize parameters
intError = [];

%% Import LDFE-SQ quadrature
dataLDFE = dlmread('LDFEQUADSALOCAL80.DAT');

% Go through each refinement level
for i = 1 : size(ref_LDFE, 2);
    
    fprintf('Local refinement level: %i \n', i);
    
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
        omegaX(j)   = dataLDFE(lineUse, 1);
        omegaY(j)   = dataLDFE(lineUse, 2);
        omegaZ(j)   = dataLDFE(lineUse, 3);
        weights(j)  = dataLDFE(lineUse, 4);
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