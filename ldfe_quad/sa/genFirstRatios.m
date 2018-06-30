% Generate LDFE-ratio quadrature
function [quadrature, rhoInitAll] = genFirstRatios(squareInfo, id)

% Retrieve square properties
midSubX      = squareInfo.midSubX;
midSubY      = squareInfo.midSubY;
minSubSubX   = squareInfo.minSubSubX;
maxSubSubX   = squareInfo.maxSubSubX;
minSubSubY   = squareInfo.minSubSubY;
maxSubSubY   = squareInfo.maxSubSubY;
surfaceArea  = squareInfo.surfaceArea;
integrations = squareInfo.integrations;
numSubSq = squareInfo.numSubSq;

% Initialize storage
xPos    = cell(numSubSq, 1);
yPos    = cell(numSubSq, 1);
gamma   = cell(numSubSq, 3);
theta   = cell(numSubSq, 3);
weights = cell(numSubSq, 1);
rhoInitAll = zeros(1, numSubSq);

% Define iteration parameters
areaEps   = 1e-15;
maxCounts = 10000;

% Go through each sub-square
for iSub = 1 : numSubSq
    
    % Upper bound guess
    ratio_a = 0.8;
    
    % Lower bound guess
    ratio_b = 0.4;
    
    % Solve initial weights for a ratio
    xPosTemp = [...
        midSubX(iSub) + (maxSubSubX{iSub}(1) - minSubSubX{iSub}(1)) * ratio_a,...
        midSubX(iSub) + (maxSubSubX{iSub}(2) - minSubSubX{iSub}(2)) * ratio_a,...
        midSubX(iSub) - (maxSubSubX{iSub}(3) - minSubSubX{iSub}(3)) * ratio_a,...
        midSubX(iSub) - (maxSubSubX{iSub}(4) - minSubSubX{iSub}(4)) * ratio_a];
    yPosTemp = [...
        midSubY(iSub) + (maxSubSubY{iSub}(1) - minSubSubY{iSub}(1)) * ratio_a,...
        midSubY(iSub) - (maxSubSubY{iSub}(2) - minSubSubY{iSub}(2)) * ratio_a,...
        midSubY(iSub) - (maxSubSubY{iSub}(3) - minSubSubY{iSub}(3)) * ratio_a,...
        midSubY(iSub) + (maxSubSubY{iSub}(4) - minSubSubY{iSub}(4)) * ratio_a];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations{iSub}, constants);
    
    % Initial residual
    res_a = (surfaceArea{iSub}(id) - weightsTemp(id)) / surfaceArea{iSub}(id);
    
    % Make sure residual is positive
    if res_a < 0
        error('ratio_a gave negative weights!');
    end
    
    % Solve new weights
    xPosTemp = [...
        midSubX(iSub) + (maxSubSubX{iSub}(1) - minSubSubX{iSub}(1)) * ratio_b,...
        midSubX(iSub) + (maxSubSubX{iSub}(2) - minSubSubX{iSub}(2)) * ratio_b,...
        midSubX(iSub) - (maxSubSubX{iSub}(3) - minSubSubX{iSub}(3)) * ratio_b,...
        midSubX(iSub) - (maxSubSubX{iSub}(4) - minSubSubX{iSub}(4)) * ratio_b];
    yPosTemp = [...
        midSubY(iSub) + (maxSubSubY{iSub}(1) - minSubSubY{iSub}(1)) * ratio_b,...
        midSubY(iSub) - (maxSubSubY{iSub}(2) - minSubSubY{iSub}(2)) * ratio_b,...
        midSubY(iSub) - (maxSubSubY{iSub}(3) - minSubSubY{iSub}(3)) * ratio_b,...
        midSubY(iSub) + (maxSubSubY{iSub}(4) - minSubSubY{iSub}(4)) * ratio_b];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations{iSub}, constants);
    
    % Solve new residual
    res_b = (surfaceArea{iSub}(id) - weightsTemp(id)) / surfaceArea{iSub}(id);
    
    % Make sure residual is negative
    if res_b > 0
        error('ratio_b is positive!');
    end
    
    % Iterate until chosen sub-square weight equals its surface area
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts                
        
        % Calculate next guess
        x = (ratio_b + ratio_a) / 2;
        
        % New weights for next guess
        xPosTemp = [midSubX(iSub) + (maxSubSubX{iSub}(1) - minSubSubX{iSub}(1)) * x,...
                    midSubX(iSub) + (maxSubSubX{iSub}(2) - minSubSubX{iSub}(2)) * x,...
                    midSubX(iSub) - (maxSubSubX{iSub}(3) - minSubSubX{iSub}(3)) * x,...
                    midSubX(iSub) - (maxSubSubX{iSub}(4) - minSubSubX{iSub}(4)) * x];
        yPosTemp = [midSubY(iSub) + (maxSubSubY{iSub}(1) - minSubSubY{iSub}(1)) * x,...
                    midSubY(iSub) - (maxSubSubY{iSub}(2) - minSubSubY{iSub}(2)) * x,...
                    midSubY(iSub) - (maxSubSubY{iSub}(3) - minSubSubY{iSub}(3)) * x,...
                    midSubY(iSub) + (maxSubSubY{iSub}(4) - minSubSubY{iSub}(4)) * x];
        [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
        constants = Basis(gammaTemp, thetaTemp);
        weightsTemp = Weights(integrations{iSub}, constants);
        
        % New residual
        res_x = (surfaceArea{iSub}(id) - weightsTemp(id)) / surfaceArea{iSub}(id);
                
        % Ensure residuals are not equal
        if res_a == res_b
            error('Error: residuals are equal');
        else            
            
            % Check for convergence
            if abs(ratio_a - ratio_b) < areaEps               
                converged = 1;            
            else
                
                % If the residual of next guess has same sign as a
                if res_x / res_a > 0
                    
                    % New a ratio
                    ratio_a = x;
                    
                    % New a residual
                    res_a = res_x;
                    
                else
                    % New b ratio
                    ratio_b = x;
                    
                    % New b residual
                    res_b = res_x;
                    
                end
                
                % Update counter
                counter = counter + 1;
                
            end
        end
    end
    
    % Check for failure
    if converged == 0
        error('Error: LDFE-ratio failed');
    end
    
    % Store quadrature data
    xPos{iSub}     = xPosTemp;
    yPos{iSub}     = yPosTemp;
    gamma{iSub, 1} = gammaTemp;
    theta{iSub, 1} = thetaTemp;
    weights{iSub}  = weightsTemp;
    rhoInitAll(iSub) = ratio_a;
    
end

% Maximum surface area error
test = zeros(numSubSq, 1);
counter = 1;
for i = 1 : numSubSq
    test(counter) = abs(weights{i}(id) - surfaceArea{i}(id)) / surfaceArea{i}(id);
    counter = counter + 1;
end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

% Rotate to other three faces
for i = 1 : numSubSq
    for j = 1 : 4
        gamma{i, 2}(j) = (pi / 2) - atan(sqrt(3) * yPos{i}(j));
        theta{i, 2}(j) = (pi / 2) - atan(sqrt(3) * xPos{i}(j) /...
            sqrt(1 + 3 * yPos{i}(j) ^ 2));
        gamma{i, 3}(j) = atan(yPos{i}(j) / xPos{i}(j));
        theta{i, 3}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos{i}(j) ^ 2 + ...
            xPos{i}(j) ^ 2)));
    end
end

% Store LDFE-center quadrature data
quadrature.gamma = gamma;
quadrature.theta = theta;
quadrature.weights = weights;

% Create and store directional cosines
counter = 1;
for i = 1 : numSubSq
    for j = 1 : 4
        for k = 1 : 3
            omegaX(counter) = cos(gamma{i, k}(j))*sin(theta{i, k}(j));
            omegaY(counter) = sin(gamma{i, k}(j))*sin(theta{i, k}(j));
            omegaZ(counter) = cos(theta{i, k}(j));
            weight(counter) = weights{i}(j);
            counter = counter + 1;
        end        
    end
end
quadrature.omegaX = omegaX;
quadrature.omegaY = omegaY;
quadrature.omegaZ = omegaZ;
quadrature.weight = weight;

end