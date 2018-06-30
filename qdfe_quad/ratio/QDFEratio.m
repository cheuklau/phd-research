% Generate QDFE-ratio quadrature
function quadrature = QDFEratio(squareInfo, iRef)

% Retrieve needed square properties
minSubX       = squareInfo.minSubX;
minSubY       = squareInfo.minSubY;
maxSubX       = squareInfo.maxSubX;
maxSubY       = squareInfo.maxSubY;
minSubSubX    = squareInfo.minSubSubX;
maxSubSubX    = squareInfo.maxSubSubX;
minSubSubY    = squareInfo.minSubSubY;
maxSubSubY    = squareInfo.maxSubSubY;
surfaceArea   = squareInfo.surfaceArea;
integrations  = squareInfo.integrations;

% Initialize storage
xPos = cell((iRef + 1) ^ 2, 1);
yPos = cell((iRef + 1) ^ 2, 1);
gamma = cell((iRef + 1) ^ 2, 3);
theta = cell((iRef + 1) ^ 2, 3);
weights = cell((iRef + 1) ^ 2, 1);

% Go through each sub-square
areaEps = 1e-12;
maxCounts = 50;
ratioFixed = 0.5; % Variable, used for non-corner sub-sub-squares
for iSub = 1 : (iRef + 1) ^ 2
    % Initial ratio guesses
    ratio = [0.6 0.3];
    ratio0 = ratio(1);
    % Solve initial weights
    subSubLengthX = (maxSubX(iSub) - minSubX(iSub)) / 3;
    subSubLengthY = (maxSubY(iSub) - minSubY(iSub)) / 3;
    xPosTemp = [maxSubSubX{iSub}(1) - subSubLengthX * ratio0, ...
        (maxSubSubX{iSub}(2) + minSubSubX{iSub}(2)) / 2, ...
        minSubSubX{iSub}(3) + subSubLengthX * ratio0, ...
        maxSubSubX{iSub}(4) - subSubLengthX * ratioFixed, ...
        (maxSubSubX{iSub}(5) + minSubSubX{iSub}(5)) / 2, ...
        minSubSubX{iSub}(6) + subSubLengthX * ratioFixed, ...
        maxSubSubX{iSub}(7) - subSubLengthX * ratio0, ...
        (maxSubSubX{iSub}(8) + minSubSubX{iSub}(8)) / 2, ...
        minSubSubX{iSub}(9) + subSubLengthX * ratio0];
    yPosTemp = [maxSubSubY{iSub}(1) - subSubLengthY * ratio0, ...
        maxSubSubY{iSub}(2) - subSubLengthY * ratioFixed, ...
        maxSubSubY{iSub}(3) - subSubLengthY * ratio0, ...
        (maxSubSubY{iSub}(4) + minSubSubY{iSub}(4)) / 2, ...
        (maxSubSubY{iSub}(5) + minSubSubY{iSub}(5)) / 2, ...
        (maxSubSubY{iSub}(6) + minSubSubY{iSub}(6)) / 2, ...
        minSubSubY{iSub}(7) + subSubLengthY * ratio0, ...
        minSubSubY{iSub}(8) + subSubLengthY * ratioFixed, ...
        minSubSubY{iSub}(9) + subSubLengthY * ratio0];
    [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
    constants = Basis(gammaTemp, thetaTemp);
    weights0 = Weights(integrations{iSub}, constants);
    % Initial difference for chosen sub-square
    res(1) = (surfaceArea{iSub}(5) - weights0(5)) / surfaceArea{iSub}(5);
    % Iterate until chosen sub-square weight equals its surface area
    counter = 0;
    converged = 0;
    while converged == 0 && counter < maxCounts
        % Next ratio guess
        ratio1 = ratio(2);
        % Solve new weights
        xPosTemp = [maxSubSubX{iSub}(1) - subSubLengthX * ratio1, ...
            (maxSubSubX{iSub}(2) + minSubSubX{iSub}(2)) / 2, ...
            minSubSubX{iSub}(3) + subSubLengthX * ratio1, ...
            maxSubSubX{iSub}(4) - subSubLengthX * ratioFixed, ...
            (maxSubSubX{iSub}(5) + minSubSubX{iSub}(5)) / 2, ...
            minSubSubX{iSub}(6) + subSubLengthX * ratioFixed, ...
            maxSubSubX{iSub}(7) - subSubLengthX * ratio1, ...
            (maxSubSubX{iSub}(8) + minSubSubX{iSub}(8)) / 2, ...
            minSubSubX{iSub}(9) + subSubLengthX* ratio1];
        yPosTemp = [maxSubSubY{iSub}(1) - subSubLengthY * ratio1, ...
            maxSubSubY{iSub}(2) - subSubLengthY * ratioFixed, ...
            maxSubSubY{iSub}(3) - subSubLengthY * ratio1, ...
            (maxSubSubY{iSub}(4) + minSubSubY{iSub}(4)) / 2, ...
            (maxSubSubY{iSub}(5) + minSubSubY{iSub}(5)) / 2, ...
            (maxSubSubY{iSub}(6) + minSubSubY{iSub}(6)) / 2, ...
            minSubSubY{iSub}(7) + subSubLengthY * ratio1, ...
            minSubSubY{iSub}(8) + subSubLengthY * ratioFixed, ...
            minSubSubY{iSub}(9) + subSubLengthY * ratio1];
        [gammaTemp, thetaTemp] = LocalToGlobal(xPosTemp, yPosTemp);
        constants = Basis(gammaTemp, thetaTemp);
        weights1 = Weights(integrations{iSub}, constants);
        res(2) = (surfaceArea{iSub}(5) - weights1(5)) / surfaceArea{iSub}(5);
        % Ensure at least one step
        if res(1) == res(2)
            ratio(2) = (ratio(1) + ratio(2)) / 2;
            ratio(1) = ratio(2);
            counter = counter + 1;
        else
            delta = (ratio(2) - ratio(1)) / (res(2) - res(1)) * res(2);
            if abs(delta) < areaEps
                ratio(2) = ratio(2) - delta;
                converged = 1;
            else
                % Update ratios
                res(1) = res(2);
                ratio(1) = ratio(2);
                ratio(2) = ratio(2) - delta;
                counter = counter + 1;
            end
        end 
    end
    if converged == 0
      %  error('get first ratio failed!');
    end
    % Store LDFE-center quadrature data
    xPos{iSub} = xPosTemp;
    yPos{iSub} = yPosTemp;
    gamma{iSub, 1} = gammaTemp;
    theta{iSub, 1} = thetaTemp;
    weights{iSub} = weights1;    
end

% Normalize weights to 4pi
tot_weight = 0;
counter = 1;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        tot_weight = tot_weight + weights{i}(j);
        counter = counter + 1;
    end
end
norm = (4 * pi / 24) / tot_weight;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        weights{i}(j) = weights{i}(j) * norm;
    end
end

% Maximum surface area error
test = zeros((iRef + 1) ^ 2, 1);
counter = 1;
for i = 1 : (iRef + 1) ^ 2
    test(counter) = abs(weights{i}(5) - surfaceArea{i}(5)) / surfaceArea{i}(5);
    counter = counter + 1;
end
fprintf('the maximum error is: %E \n', max(test));
fprintf('the average error is: %E \n', mean(test));

% Rotate to other three faces
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9        
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

end