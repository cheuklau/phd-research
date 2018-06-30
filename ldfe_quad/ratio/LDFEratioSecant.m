% Generate LDFE-ratio quadrature
function quadrature = LDFEratioSecant(squareInfo, iRef, id)

% Retrieve square properties
maxSubX      = squareInfo.maxSubX;
maxSubY      = squareInfo.maxSubY;
midSubX      = squareInfo.midSubX;
midSubY      = squareInfo.midSubY;
surfaceArea  = squareInfo.surfaceArea;
integrations = squareInfo.integrations;

% Initialize storage
xPos    = cell((iRef + 1) ^ 2, 1);
yPos    = cell((iRef + 1) ^ 2, 1);
gamma   = cell((iRef + 1) ^ 2, 3);
theta   = cell((iRef + 1) ^ 2, 3);
weights = cell((iRef + 1) ^ 2, 1);

% Define iteration parameters
areaEps   = 1e-11;
maxCounts = 1000;

% Go through each sub-square
for iSub = 1 : (iRef + 1) ^ 2
    disp(iSub);
    % Initial ratio guesses
    ratio0 = 1/sqrt(3);
    %ratio0 = 0.6;
    
    % Next ratio guess
    %ratio1 = 0.4;
    ratio1 = 0.57;
    
    % Plot residual as a function of ratio
    
    ratio_test = linspace(0, 1, 100000);
    counter = 1;
    for i = 1 : size(ratio_test, 2)
        % Solve weight of current residual
        subSubLengthX = maxSubX(iSub) - midSubX(iSub);
        subSubLengthY = maxSubY(iSub) - midSubY(iSub);
        xPosTemp = [midSubX(iSub) + subSubLengthX * ratio_test(i),...
            midSubX(iSub) + subSubLengthX * ratio_test(i),...
            midSubX(iSub) - subSubLengthX * ratio_test(i),...
            midSubX(iSub) - subSubLengthX * ratio_test(i)];
        yPosTemp = [midSubY(iSub) + subSubLengthY * ratio_test(i),...
            midSubY(iSub) - subSubLengthY * ratio_test(i),...
            midSubY(iSub) - subSubLengthY * ratio_test(i),...
            midSubY(iSub) + subSubLengthY * ratio_test(i)];
        [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
        constants = Basis(gammaTemp, thetaTemp);
        weightsTemp = Weights(integrations{iSub}, constants);
        
        % Initial residual        
        res_temp = (surfaceArea{iSub}(id) - weightsTemp(id)) / surfaceArea{iSub}(id);
        if res_temp < 0
            res_plot(counter) = abs(res_temp);
            interval_plot(counter) = ratio_test(i);
            counter = counter + 1;
        end
        
    end    

    % Solve initial weights
    subSubLengthX = maxSubX(iSub) - midSubX(iSub);
    subSubLengthY = maxSubY(iSub) - midSubY(iSub);
    xPosTemp = [midSubX(iSub) + subSubLengthX * ratio0,...
        midSubX(iSub) + subSubLengthX * ratio0,...
        midSubX(iSub) - subSubLengthX * ratio0,...
        midSubX(iSub) - subSubLengthX * ratio0];
    yPosTemp = [midSubY(iSub) + subSubLengthY * ratio0,...
        midSubY(iSub) - subSubLengthY * ratio0,...
        midSubY(iSub) - subSubLengthY * ratio0,...
        midSubY(iSub) + subSubLengthY * ratio0];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weightsTemp = Weights(integrations{iSub}, constants);
    
    % Initial residual
    res0 = (surfaceArea{iSub}(id) - weightsTemp(id)) / surfaceArea{iSub}(id);
    
    % Iterate until chosen sub-square weight equals its surface area
    counter   = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts
        
        % Solve new weights
        xPosTemp = [midSubX(iSub) + subSubLengthX * ratio1,...
            midSubX(iSub) + subSubLengthX * ratio1,...
            midSubX(iSub) - subSubLengthX * ratio1,...
            midSubX(iSub) - subSubLengthX * ratio1];
        yPosTemp = [midSubY(iSub) + subSubLengthY * ratio1,...
            midSubY(iSub) - subSubLengthY * ratio1,...
            midSubY(iSub) - subSubLengthY * ratio1,...
            midSubY(iSub) + subSubLengthY * ratio1];
        [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
        constants = Basis(gammaTemp, thetaTemp);
        weightsTemp = Weights(integrations{iSub}, constants);
        
        % Solve new residual
        res1 = (surfaceArea{iSub}(id) - weightsTemp(id)) / surfaceArea{iSub}(id);
        
        % Ensure residuals are not equal
        if res0 == res1
            error('Error: residuals are equal');
        else
            
            % Calculate delta and check for convergence
            delta = (ratio1 - ratio0) / (res1 - res0) * res1;
            if abs(delta) < areaEps
                ratio1 = ratio1 - delta;
                converged = 1;
                ratio_store(iSub) = ratio1;
                res_store(iSub) = abs(res1);
            else
                
                % Update ratios
                res0 = res1;
                ratio0 = ratio1;
                ratio1 = ratio1 - delta;
                counter = counter + 1;
            end
        end
    end
    if converged == 0
        error('Error: LDFE-ratio failed');
    end
    
    % Store quadrature data
    xPos{iSub}     = xPosTemp;
    yPos{iSub}     = yPosTemp;
    gamma{iSub, 1} = gammaTemp;
    theta{iSub, 1} = thetaTemp;
    weights{iSub}  = weightsTemp;
    
end

% Maximum surface area error
test = zeros((iRef + 1) ^ 2, 1);
counter = 1;
for i = 1 : (iRef + 1) ^ 2
    test(counter) = abs(weights{i}(id) - surfaceArea{i}(id)) / surfaceArea{i}(id);
    counter = counter + 1;
end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

% Rotate to other three faces
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 4
        gamma{i, 2}(j) = (pi / 2) - atan(sqrt(3) * yPos{i}(j));
        theta{i, 2}(j) = (pi / 2) - atan(sqrt(3) * xPos{i}(j) / sqrt(1 + 3 * yPos{i}(j) ^ 2));
        gamma{i, 3}(j) = atan(yPos{i}(j) / xPos{i}(j));
        theta{i, 3}(j) = (pi / 2) - atan(1 / (sqrt(3) * sqrt(yPos{i}(j) ^ 2 + xPos{i}(j) ^ 2)));
    end
end

% Store LDFE-center quadrature data
quadrature.gamma   = gamma;
quadrature.theta   = theta;
quadrature.weights = weights;

end