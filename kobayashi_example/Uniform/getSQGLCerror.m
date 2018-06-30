% Master file for cube-based QDFE quadrature sets
function [intErrorSQGLC, meshSize] = getSQGLCerror(v_top, v_bot, u_tip, flip, refLevels)

% Number of lines
for i = 1 : size(refLevels, 2)
    
    numLines(i) = refLevels(i) ^ 2;
    
end

% Initialize relative error
intErrorSQGLC = [];

% Import GLC data
dataSQGLC = dlmread('SQGLC.DAT');

% Go through each refinement level
for i = 1 : size(refLevels, 2);
    
    fprintf('Currently on refinement level: %i \n', refLevels(i));
    
    % Starting line
    lineStart = 1;
    for j = 1 : i - 1
        lineStart = lineStart + numLines(j);
    end
    
    % Store data values    
    omegaX = zeros(1, numLines(i));
    omegaY = zeros(1, numLines(i));
    omegaZ = zeros(1, numLines(i));
    weight = zeros(1, numLines(i));
    weightTotal = 0;
    for j = 1 : numLines(i)
        lineUse = lineStart + j - 1;
        omegaX(j) = dataSQGLC(lineUse, 1);
        omegaY(j) = dataSQGLC(lineUse, 2);
        omegaZ(j) = dataSQGLC(lineUse, 3);
        weight(j) = dataSQGLC(lineUse, 4);
        weightTotal = weightTotal + weight(j);
    end
    
    % Normalize weights
    for j = 1 : numLines(i)
       weight(j) = weight(j) * (4 * pi / 8) / weightTotal; 
    end
    
    % Calculate error
    intErrorSQGLC = SphInt(omegaX, omegaY, omegaZ, weight, intErrorSQGLC, v_top, v_bot, u_tip, flip);
    
    % Angular mesh size
    meshSize(i) = 1 / sqrt(8 * numLines(i));
        
end

end