% Calculate basic properties of the square with
% equal arc-length sub-squares and sub-sub-squares
function squareInfo = SquareProp(ref0, ref1, perInc)

%% Cone of angle

% Gamma limits
gammaMinCone = 1.17663374;
gammaMaxCone = 1.31657578;
gammaDif = perInc * (gammaMaxCone - gammaMinCone);
gammaMinCone = gammaMinCone - gammaDif;
gammaMaxCone = gammaMaxCone + gammaDif;

% Theta limits
thetaMinCone = 1.19093456;
thetaMaxCone = 1.33176512;
thetaDif = perInc * (thetaMaxCone - thetaMinCone);
thetaMinCone = thetaMinCone - thetaDif;
thetaMaxCone = thetaMaxCone + thetaDif;

%% Generate uniform quadratures for initial refinement level

% Initialize storage
minSubX = cell(ref1 - ref0 + 1, 1);
maxSubX = cell(ref1 - ref0 + 1, 1);
minSubY = cell(ref1 - ref0 + 1, 1);
maxSubY = cell(ref1 - ref0 + 1, 1);
midSubX = cell(ref1 - ref0 + 1, 1);
midSubY = cell(ref1 - ref0 + 1, 1);
minSubSubX = cell(ref1 - ref0 + 1, 1);
maxSubSubX = cell(ref1 - ref0 + 1, 1);
minSubSubY = cell(ref1 - ref0 + 1, 1);
maxSubSubY = cell(ref1 - ref0 + 1, 1);
midSubSubX = cell(ref1 - ref0 + 1, 1);
midSubSubY = cell(ref1 - ref0 + 1, 1);
xPos   = cell(ref1 - ref0 + 1, 1);
yPos   = cell(ref1 - ref0 + 1, 1);
numSubSq = zeros(ref1 - ref0 + 1, 1);

% Go through each refinement level
for iRef = ref0 : ref1
    
    % Direct refinement equivalent
    iRefDir = sqrt(2 ^ (iRef * 2)) - 1;
    [minSubX{iRef}, maxSubX{iRef}, ...
     minSubY{iRef}, maxSubY{iRef}, ...
     midSubX{iRef}, midSubY{iRef}, ...
     midSubSubX{iRef}, midSubSubY{iRef}, ...
     minSubSubX{iRef}, maxSubSubX{iRef}, ...
     minSubSubY{iRef}, maxSubSubY{iRef}, ...
     xPos{iRef}, yPos{iRef}, ...
     numSubSq(iRef)] = calcUniform(iRefDir);    
end

%% Generate locally-refined quadrature for first face

% Initialize storage
minSubX1    = minSubX{1};
maxSubX1    = maxSubX{1};
minSubY1    = minSubY{1};
maxSubY1    = maxSubY{1};
midSubX1    = midSubX{1};
midSubY1    = midSubY{1};
minSubSubX1 = minSubSubX{1};
maxSubSubX1 = maxSubSubX{1};
minSubSubY1 = minSubSubY{1};
maxSubSubY1 = maxSubSubY{1};
midSubSubX1 = midSubSubX{1};
midSubSubY1 = midSubSubY{1};
xPos1       = xPos{1};
yPos1       = yPos{1};
numSubSq1   = numSubSq(1);

% Go through each refinement level
for i = ref0 : ref1 - 1
   
    % Temporary storage for number of sub-squares
    numSubSq1tmp = numSubSq1;
    
    % Go through each sub-square    
    for j = 1 : numSubSq1
        
        gamma = [];
        theta = [];
        
        % Check if sub-square is in cone of angle
        gamma(1) = atan(sqrt(3) * minSubX1(j));
        gamma(2) = atan(sqrt(3) * maxSubX1(j));
        gammaMin = min(gamma);
        gammaMax = max(gamma);
        theta(1) = pi / 2 - atan(sqrt(3 / (1 + 3 * minSubX1(j) ^ 2)) * maxSubY1(j));
        theta(2) = pi / 2 - atan(sqrt(3 / (1 + 3 * minSubX1(j) ^ 2)) * minSubY1(j));
        theta(3) = pi / 2 - atan(sqrt(3 / (1 + 3 * maxSubX1(j) ^ 2)) * maxSubY1(j));
        theta(4) = pi / 2 - atan(sqrt(3 / (1 + 3 * maxSubX1(j) ^ 2)) * minSubY1(j));
        thetaMin = min(theta);
        thetaMax = max(theta);
        
        if ... corner case
                (((gammaMin >= gammaMinCone && gammaMin <= gammaMaxCone)   || ...
                  (gammaMax >= gammaMinCone && gammaMax <= gammaMaxCone))  && ...
                 ((thetaMin >= thetaMinCone && thetaMin <= thetaMaxCone)   || ...
                  (thetaMax >= thetaMinCone && thetaMax <= thetaMaxCone))) || ...
                ... side case along theta
                (((gammaMin >= gammaMinCone && gammaMin <= gammaMaxCone)   || ...
                  (gammaMax >= gammaMinCone && gammaMax <= gammaMaxCone))  && ...
                  (thetaMin <= thetaMinCone && thetaMax >= thetaMaxCone))  || ...
                ... side case along gamma
                (((thetaMin >= thetaMinCone && thetaMin <= thetaMaxCone)   || ...
                  (thetaMax >= thetaMinCone && thetaMax <= thetaMaxCone))  && ...
                  (gammaMin <= gammaMinCone && gammaMax >= gammaMaxCone))  || ...
                ... center case
                (gammaMin <= gammaMinCone && gammaMax >= gammaMaxCone && ...
                 thetaMin <= thetaMinCone && thetaMax >= thetaMaxCone)
            
            % Upper-right refined sub-square
            xPosTmp = (xPos1(j) - 1) * 2 + 2; % x-index of refined sub-square
            yPosTmp = (yPos1(j) - 1) * 2 + 2; % y-index of refined sub-square
            index = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX1(j) = minSubX{i + 1}(index);
            maxSubX1(j) = maxSubX{i + 1}(index);
            minSubY1(j) = minSubY{i + 1}(index);
            maxSubY1(j) = maxSubY{i + 1}(index);
            midSubX1(j) = midSubX{i + 1}(index);
            midSubY1(j) = midSubY{i + 1}(index);
            midSubSubX1{j} = midSubSubX{i + 1}{index};
            midSubSubY1{j} = midSubSubY{i + 1}{index};
            minSubSubX1{j} = minSubSubX{i + 1}{index};
            maxSubSubX1{j} = maxSubSubX{i + 1}{index};
            minSubSubY1{j} = minSubSubY{i + 1}{index};
            maxSubSubY1{j} = maxSubSubY{i + 1}{index};            
            
            % Upper-left refined sub-square
            xPosTmp = (xPos1(j) - 1) * 2 + 1;
            yPosTmp = (yPos1(j) - 1) * 2 + 2;
            index  = numSubSq1tmp + 1;
            index2 = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX1(index) = minSubX{i + 1}(index2);
            maxSubX1(index) = maxSubX{i + 1}(index2);
            minSubY1(index) = minSubY{i + 1}(index2);
            maxSubY1(index) = maxSubY{i + 1}(index2);
            midSubX1(index) = midSubX{i + 1}(index2);
            midSubY1(index) = midSubY{i + 1}(index2);
            midSubSubX1{index} = midSubSubX{i + 1}{index2};
            midSubSubY1{index} = midSubSubY{i + 1}{index2};
            minSubSubX1{index} = minSubSubX{i + 1}{index2};
            maxSubSubX1{index} = maxSubSubX{i + 1}{index2};
            minSubSubY1{index} = minSubSubY{i + 1}{index2};
            maxSubSubY1{index} = maxSubSubY{i + 1}{index2};
            xPos1(index) = xPosTmp;
            yPos1(index) = yPosTmp;
            
            % Lower-right refined sub-square
            xPosTmp = (xPos1(j) - 1) * 2 + 2;
            yPosTmp = (yPos1(j) - 1) * 2 + 1;
            index = numSubSq1tmp + 2;
            index2 = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX1(index) = minSubX{i + 1}(index2);
            maxSubX1(index) = maxSubX{i + 1}(index2);
            minSubY1(index) = minSubY{i + 1}(index2);
            maxSubY1(index) = maxSubY{i + 1}(index2);
            midSubX1(index) = midSubX{i + 1}(index2);
            midSubY1(index) = midSubY{i + 1}(index2);
            midSubSubX1{index} = midSubSubX{i + 1}{index2};
            midSubSubY1{index} = midSubSubY{i + 1}{index2};
            minSubSubX1{index} = minSubSubX{i + 1}{index2};
            maxSubSubX1{index} = maxSubSubX{i + 1}{index2};
            minSubSubY1{index} = minSubSubY{i + 1}{index2};
            maxSubSubY1{index} = maxSubSubY{i + 1}{index2};
            xPos1(index) = xPosTmp;
            yPos1(index) = yPosTmp;
            
            % Lower-left refined sub-square
            xPosTmp = (xPos1(j) - 1) * 2 + 1;
            yPosTmp = (yPos1(j) - 1) * 2 + 1;
            index = numSubSq1tmpt + 3;
            index2 = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX1(index) = minSubX{i + 1}(index2);
            maxSubX1(index) = maxSubX{i + 1}(index2);
            minSubY1(index) = minSubY{i + 1}(index2);
            maxSubY1(index) = maxSubY{i + 1}(index2);
            midSubX1(index) = midSubX{i + 1}(index2);
            midSubY1(index) = midSubY{i + 1}(index2);
            midSubSubX1{index} = midSubSubX{i + 1}{index2};
            midSubSubY1{index} = midSubSubY{i + 1}{index2};
            minSubSubX1{index} = minSubSubX{i + 1}{index2};
            maxSubSubX1{index} = maxSubSubX{i + 1}{index2};
            minSubSubY1{index} = minSubSubY{i + 1}{index2};
            maxSubSubY1{index} = maxSubSubY{i + 1}{index2};
            xPos1(index) = xPosTmp;
            yPos1(index) = yPosTmp;
            
            % Update x-index and y-index of current sub-square
            xPos1(j) = (xPos1(j) - 1) * 2 + 2;
            yPos1(j) = (yPos1(j) - 1) * 2 + 2;
            
            % Update number of sub-squares in temporary storage
            numSubSq1tmp = numSubSq1tmp + 3;
            
        end
                                
    end
    
    % Update number of sub-squares
    numSubSq1 = numSubSq1tmp;
    
end

%% Generate locally-refined quadrature for second face

% Initialize storage
minSubX2    = minSubX{1};
maxSubX2    = maxSubX{1};
minSubY2    = minSubY{1};
maxSubY2    = maxSubY{1};
midSubX2    = midSubX{1};
midSubY2    = midSubY{1};
minSubSubX2 = minSubSubX{1};
maxSubSubX2 = maxSubSubX{1};
minSubSubY2 = minSubSubY{1};
maxSubSubY2 = maxSubSubY{1};
midSubSubX2 = midSubSubX{1};
midSubSubY2 = midSubSubY{1};
xPos2       = xPos{1};
yPos2       = yPos{1};
numSubSq2   = numSubSq(1);

% Go through each refinement level
for i = ref0 : ref1 - 1
   
    % Direct refinement index 
    refDir = sqrt(2 ^ (i * 2)) - 1;
    
    % Temporary storage for number of sub-squares
    numSubSq2tmp = numSubSq2;
    
    % Go through each sub-square    
    for j = 1 : numSubSq2
        
        gamma = [];
        theta = [];
        
        % Check if sub-square is in cone of angle
        gamma(1)  = pi / 2 - atan(sqrt(3) * maxSubY2(j));
        gamma(2)  = pi / 2 - atan(sqrt(3) * minSubY2(j));
        gammaMin = min(gamma);
        gammaMax = max(gamma);
        theta(1) = pi / 2 - atan(sqrt(3 / (1 + 3 * minSubY2(j) ^ 2)) * maxSubX2(j));
        theta(2) = pi / 2 - atan(sqrt(3 / (1 + 3 * minSubY2(j) ^ 2)) * minSubX2(j));
        theta(3) = pi / 2 - atan(sqrt(3 / (1 + 3 * maxSubY2(j) ^ 2)) * maxSubX2(j));
        theta(4) = pi / 2 - atan(sqrt(3 / (1 + 3 * maxSubY2(j) ^ 2)) * minSubX2(j));
        thetaMin = min(theta);
        thetaMax = max(theta);
        
        if ... corner case
                (((gammaMin >= gammaMinCone && gammaMin <= gammaMaxCone)   || ...
                  (gammaMax >= gammaMinCone && gammaMax <= gammaMaxCone))  && ...
                 ((thetaMin >= thetaMinCone && thetaMin <= thetaMaxCone)   || ...
                  (thetaMax >= thetaMinCone && thetaMax <= thetaMaxCone))) || ...
                ... side case along theta
                (((gammaMin >= gammaMinCone && gammaMin <= gammaMaxCone)   || ...
                  (gammaMax >= gammaMinCone && gammaMax <= gammaMaxCone))  && ...
                  (thetaMin <= thetaMinCone && thetaMax >= thetaMaxCone))  || ...
                ... side case along gamma
                (((thetaMin >= thetaMinCone && thetaMin <= thetaMaxCone)   || ...
                  (thetaMax >= thetaMinCone && thetaMax <= thetaMaxCone))  && ...
                  (gammaMin <= gammaMinCone && gammaMax >= gammaMaxCone))  || ...
                ... center case
                  (gammaMin <= gammaMinCone && gammaMax >= gammaMaxCone && ...
                   thetaMin <= thetaMinCone && thetaMax >= thetaMaxCone)
            
            % Upper-right refined sub-square
            xPosTmp = (xPos2(j) - 1) * 2 + 2; % x-index of refined sub-square
            yPosTmp = (yPos2(j) - 1) * 2 + 2; % y-index of refined sub-square
            index = (yPosTmp - 1) * (refDir + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX2(j) = minSubX{i + 1}(index);
            maxSubX2(j) = maxSubX{i + 1}(index);
            minSubY2(j) = minSubY{i + 1}(index);
            maxSubY2(j) = maxSubY{i + 1}(index);
            midSubX2(j) = midSubX{i + 1}(index);
            midSubY2(j) = midSubY{i + 1}(index);
            midSubSubX2{j} = midSubSubX{i + 1}{index};
            midSubSubY2{j} = midSubSubY{i + 1}{index};
            minSubSubX2{j} = minSubSubX{i + 1}{index};
            maxSubSubX2{j} = maxSubSubX{i + 1}{index};
            minSubSubY2{j} = minSubSubY{i + 1}{index};
            maxSubSubY2{j} = maxSubSubY{i + 1}{index};            
            
            % Upper-left refined sub-square
            xPosTmp = (xPos2(j) - 1) * 2 + 1;
            yPosTmp = (yPos2(j) - 1) * 2 + 2;
            index  = numSubSq2tmp + 1;
            index2 = (yPosTmp - 1) * (refDir + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX2(index) = minSubX{i + 1}(index2);
            maxSubX2(index) = maxSubX{i + 1}(index2);
            minSubY2(index) = minSubY{i + 1}(index2);
            maxSubY2(index) = maxSubY{i + 1}(index2);
            midSubX2(index) = midSubX{i + 1}(index2);
            midSubY2(index) = midSubY{i + 1}(index2);
            midSubSubX2{index} = midSubSubX{i + 1}{index2};
            midSubSubY2{index} = midSubSubY{i + 1}{index2};
            minSubSubX2{index} = minSubSubX{i + 1}{index2};
            maxSubSubX2{index} = maxSubSubX{i + 1}{index2};
            minSubSubY2{index} = minSubSubY{i + 1}{index2};
            maxSubSubY2{index} = maxSubSubY{i + 1}{index2};
            xPos2(index) = xPosTmp;
            yPos2(index) = yPosTmp;
            
            % Lower-right refined sub-square
            xPosTmp = (xPos2(j) - 1) * 2 + 2;
            yPosTmp = (yPos2(j) - 1) * 2 + 1;
            index = numSubSq2tmp + 2;
            index2 = (yPosTmp - 1) * (refDir + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX2(index) = minSubX{i + 1}(index2);
            maxSubX2(index) = maxSubX{i + 1}(index2);
            minSubY2(index) = minSubY{i + 1}(index2);
            maxSubY2(index) = maxSubY{i + 1}(index2);
            midSubX2(index) = midSubX{i + 1}(index2);
            midSubY2(index) = midSubY{i + 1}(index2);
            midSubSubX2{index} = midSubSubX{i + 1}{index2};
            midSubSubY2{index} = midSubSubY{i + 1}{index2};
            minSubSubX2{index} = minSubSubX{i + 1}{index2};
            maxSubSubX2{index} = maxSubSubX{i + 1}{index2};
            minSubSubY2{index} = minSubSubY{i + 1}{index2};
            maxSubSubY2{index} = maxSubSubY{i + 1}{index2};
            xPos2(index) = xPosTmp;
            yPos2(index) = yPosTmp;
            
            % Lower-left refined sub-square
            xPosTmp = (xPos2(j) - 1) * 2 + 1;
            yPosTmp = (yPos2(j) - 1) * 2 + 1;
            index = numSubSq2tmp + 3;
            index2 = (yPosTmp - 1) * (refDir + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX2(index) = minSubX{i + 1}(index2);
            maxSubX2(index) = maxSubX{i + 1}(index2);
            minSubY2(index) = minSubY{i + 1}(index2);
            maxSubY2(index) = maxSubY{i + 1}(index2);
            midSubX2(index) = midSubX{i + 1}(index2);
            midSubY2(index) = midSubY{i + 1}(index2);
            midSubSubX2{index} = midSubSubX{i + 1}{index2};
            midSubSubY2{index} = midSubSubY{i + 1}{index2};
            minSubSubX2{index} = minSubSubX{i + 1}{index2};
            maxSubSubX2{index} = maxSubSubX{i + 1}{index2};
            minSubSubY2{index} = minSubSubY{i + 1}{index2};
            maxSubSubY2{index} = maxSubSubY{i + 1}{index2};
            xPos2(index) = xPosTmp;
            yPos2(index) = yPosTmp;
            
            % Update x-index and y-index of current sub-square
            xPos2(j) = (xPos2(j) - 1) * 2 + 2;
            yPos2(j) = (yPos2(j) - 1) * 2 + 2;
            
            % Update number of sub-squares in temporary storage
            numSubSq2tmp = numSubSq2tmp + 3;
            
        end
                                
    end
    
    % Update number of sub-squares
    numSubSq2 = numSubSq2tmp;
    
end

%% Generate locally-refined quadrature for third face

% Initialize storage
minSubX3    = minSubX{1};
maxSubX3    = maxSubX{1};
minSubY3    = minSubY{1};
maxSubY3    = maxSubY{1};
midSubX3    = midSubX{1};
midSubY3    = midSubY{1};
minSubSubX3 = minSubSubX{1};
maxSubSubX3 = maxSubSubX{1};
minSubSubY3 = minSubSubY{1};
maxSubSubY3 = maxSubSubY{1};
midSubSubX3 = midSubSubX{1};
midSubSubY3 = midSubSubY{1};
xPos3       = xPos{1};
yPos3       = yPos{1};
numSubSq3   = numSubSq(1);

% Go through each refinement level
for i = ref0 : ref1 - 1
   
    % Temporary storage for number of sub-squares
    numSubSq3tmp = numSubSq3;
    
    % Go through each sub-square    
    for j = 1 : numSubSq3
        
        gamma = [];
        theta = [];
        
        % Check if sub-square is in cone of angle
        gamma(1) = atan(minSubY3(j) / maxSubX3(j));
        gamma(2) = atan(minSubY3(j) / minSubX3(j));
        gamma(3) = atan(maxSubY3(j) / maxSubX3(j));
        gamma(4) = atan(maxSubY3(j) / minSubX3(j));
        gammaMin = min(gamma);
        gammaMax = max(gamma);
        theta(1) = pi / 2 - atan(1 / sqrt(3 * (minSubY3(j) ^ 2 + maxSubX3(j) ^ 2)));
        theta(2) = pi / 2 - atan(1 / sqrt(3 * (minSubY3(j) ^ 2 + minSubX3(j) ^ 2)));
        theta(3) = pi / 2 - atan(1 / sqrt(3 * (maxSubY3(j) ^ 2 + maxSubX3(j) ^ 2)));
        theta(4) = pi / 2 - atan(1 / sqrt(3 * (maxSubY3(j) ^ 2 + minSubX3(j) ^ 2)));
        thetaMin = min(theta);
        thetaMax = max(theta);
        
        if ... corner case
                (((gammaMin >= gammaMinCone && gammaMin <= gammaMaxCone)   || ...
                  (gammaMax >= gammaMinCone && gammaMax <= gammaMaxCone))  && ...
                 ((thetaMin >= thetaMinCone && thetaMin <= thetaMaxCone)   || ...
                  (thetaMax >= thetaMinCone && thetaMax <= thetaMaxCone))) || ...
                ... side case along theta
                (((gammaMin >= gammaMinCone && gammaMin <= gammaMaxCone)   || ...
                  (gammaMax >= gammaMinCone && gammaMax <= gammaMaxCone))  && ...
                  (thetaMin <= thetaMinCone && thetaMax >= thetaMaxCone))  || ...
                ... side case along gamma
                (((thetaMin >= thetaMinCone && thetaMin <= thetaMaxCone)   || ...
                  (thetaMax >= thetaMinCone && thetaMax <= thetaMaxCone))  && ...
                  (gammaMin <= gammaMinCone && gammaMax >= gammaMaxCone))  || ...
                ... center case
                  (gammaMin <= gammaMinCone && gammaMax >= gammaMaxCone && ...
                   thetaMin <= thetaMinCone && thetaMax >= thetaMaxCone)
            
            % Upper-right refined sub-square
            xPosTmp = (xPos3(j) - 1) * 2 + 2; % x-index of refined sub-square
            yPosTmp = (yPos3(j) - 1) * 2 + 2; % y-index of refined sub-square
            index = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX3(j) = minSubX{i + 1}(index);
            maxSubX3(j) = maxSubX{i + 1}(index);
            minSubY3(j) = minSubY{i + 1}(index);
            maxSubY3(j) = maxSubY{i + 1}(index);
            midSubX3(j) = midSubX{i + 1}(index);
            midSubY3(j) = midSubY{i + 1}(index);
            midSubSubX3{j} = midSubSubX{i + 1}{index};
            midSubSubY3{j} = midSubSubY{i + 1}{index};
            minSubSubX3{j} = minSubSubX{i + 1}{index};
            maxSubSubX3{j} = maxSubSubX{i + 1}{index};
            minSubSubY3{j} = minSubSubY{i + 1}{index};
            maxSubSubY3{j} = maxSubSubY{i + 1}{index};            
            
            % Upper-left refined sub-square
            xPosTmp = (xPos3(j) - 1) * 2 + 1;
            yPosTmp = (yPos3(j) - 1) * 2 + 2;
            index  = numSubSq3tmp + 1;
            index2 = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX3(index) = minSubX{i + 1}(index2);
            maxSubX3(index) = maxSubX{i + 1}(index2);
            minSubY3(index) = minSubY{i + 1}(index2);
            maxSubY3(index) = maxSubY{i + 1}(index2);
            midSubX3(index) = midSubX{i + 1}(index2);
            midSubY3(index) = midSubY{i + 1}(index2);
            midSubSubX3{index} = midSubSubX{i + 1}{index2};
            midSubSubY3{index} = midSubSubY{i + 1}{index2};
            minSubSubX3{index} = minSubSubX{i + 1}{index2};
            maxSubSubX3{index} = maxSubSubX{i + 1}{index2};
            minSubSubY3{index} = minSubSubY{i + 1}{index2};
            maxSubSubY3{index} = maxSubSubY{i + 1}{index2};
            xPos3(index) = xPosTmp;
            yPos3(index) = yPosTmp;
            
            % Lower-right refined sub-square
            xPosTmp = (xPos3(j) - 1) * 2 + 2;
            yPosTmp = (yPos3(j) - 1) * 2 + 1;
            index = numSubSq3tmp + 2;
            index2 = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX3(index) = minSubX{i + 1}(index2);
            maxSubX3(index) = maxSubX{i + 1}(index2);
            minSubY3(index) = minSubY{i + 1}(index2);
            maxSubY3(index) = maxSubY{i + 1}(index2);
            midSubX3(index) = midSubX{i + 1}(index2);
            midSubY3(index) = midSubY{i + 1}(index2);
            midSubSubX3{index} = midSubSubX{i + 1}{index2};
            midSubSubY3{index} = midSubSubY{i + 1}{index2};
            minSubSubX3{index} = minSubSubX{i + 1}{index2};
            maxSubSubX3{index} = maxSubSubX{i + 1}{index2};
            minSubSubY3{index} = minSubSubY{i + 1}{index2};
            maxSubSubY3{index} = maxSubSubY{i + 1}{index2};
            xPos3(index) = xPosTmp;
            yPos3(index) = yPosTmp;
            
            % Lower-left refined sub-square
            xPosTmp = (xPos3(j) - 1) * 2 + 1;
            yPosTmp = (yPos3(j) - 1) * 2 + 1;
            index = numSubSq3tmp + 3;
            index2 = (yPosTmp - 1) * (i + 1) * 2 + (xPosTmp - 1) + 1;
            minSubX3(index) = minSubX{i + 1}(index2);
            maxSubX3(index) = maxSubX{i + 1}(index2);
            minSubY3(index) = minSubY{i + 1}(index2);
            maxSubY3(index) = maxSubY{i + 1}(index2);
            midSubX3(index) = midSubX{i + 1}(index2);
            midSubY3(index) = midSubY{i + 1}(index2);
            midSubSubX3{index} = midSubSubX{i + 1}{index2};
            midSubSubY3{index} = midSubSubY{i + 1}{index2};
            minSubSubX3{index} = minSubSubX{i + 1}{index2};
            maxSubSubX3{index} = maxSubSubX{i + 1}{index2};
            minSubSubY3{index} = minSubSubY{i + 1}{index2};
            maxSubSubY3{index} = maxSubSubY{i + 1}{index2};
            xPos3(index) = xPosTmp;
            yPos3(index) = yPosTmp;
            
            % Update x-index and y-index of current sub-square
            xPos3(j) = (xPos3(j) - 1) * 2 + 2;
            yPos3(j) = (yPos3(j) - 1) * 2 + 2;
            
            % Update number of sub-squares in temporary storage
            numSubSq3tmp = numSubSq3tmp + 3;
            
        end
                                
    end
    
    % Update number of sub-squares
    numSubSq3 = numSubSq3tmp;
    
end

%% Sub-sub-square surface areas

% Gauss quadrature order
N = 16;

% Keeps track of total weight
totWeight = 0;

% Weights for face 1 (y-z plane)
surfaceArea1  = cell(numSubSq1,  1);
for k = 1 : numSubSq1
    surfaceArea1{k} = zeros(4, 1);
    for m = 1 : 4
        [fx, wx] = lgwt(N, minSubSubX1{k}(m), maxSubSubX1{k}(m));
        [fy, wy] = lgwt(N, minSubSubY1{k}(m), maxSubSubY1{k}(m));
        for i = 1 : N
            for j = 1 : N
                surfaceArea1{k}(m) = surfaceArea1{k}(m) + ...
                    (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2))) * wx(i) * wy(j);
            end
        end
        totWeight = totWeight + surfaceArea1{k}(m);
    end
end

% Weights for face 2 (x-z plane)
surfaceArea2  = cell(numSubSq2,  1);
for k = 1 : numSubSq2
    surfaceArea2{k} = zeros(4, 1);
    for m = 1 : 4
        [fx, wx] = lgwt(N, minSubSubX2{k}(m), maxSubSubX2{k}(m));
        [fy, wy] = lgwt(N, minSubSubY2{k}(m), maxSubSubY2{k}(m));
        for i = 1 : N
            for j = 1 : N
                surfaceArea2{k}(m) = surfaceArea2{k}(m) + ...
                    (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2))) * wx(i) * wy(j);
            end
        end
        totWeight = totWeight + surfaceArea2{k}(m);
    end    
end

% Weights for face 3 (x-y plane)
surfaceArea3  = cell(numSubSq3,  1);
for k = 1 : numSubSq3
    surfaceArea3{k} = zeros(4, 1);
    for m = 1 : 4
        [fx, wx] = lgwt(N, minSubSubX3{k}(m), maxSubSubX3{k}(m));
        [fy, wy] = lgwt(N, minSubSubY3{k}(m), maxSubSubY3{k}(m));
        for i = 1 : N
            for j = 1 : N
                surfaceArea3{k}(m) = surfaceArea3{k}(m) + ...
                    (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2))) * wx(i) * wy(j);
            end
        end
        totWeight = totWeight + surfaceArea3{k}(m);
    end
end

% Display total weight error vs 4pi
fprintf('the total error for the SA: %E \n', abs(4 * pi - totWeight * 8) / (4 * pi));

%% Basis function integration over sub-square

% Integrations for face 1 (y-z plane)
integrations1 = cell(numSubSq1, 1);
for counter = 1 : numSubSq1
    xRange = [minSubX1(counter) maxSubX1(counter)];
    yRange = [minSubY1(counter) maxSubY1(counter)];
    integrations1{counter} = Integrations(xRange, yRange);
end

% Integrations for face 2 (x-z plane)
integrations2 = cell(numSubSq2, 1);
for counter = 1 : numSubSq2
    xRange = [minSubX2(counter) maxSubX2(counter)];
    yRange = [minSubY2(counter) maxSubY2(counter)];
    integrations2{counter} = Integrations(xRange, yRange);
end

% Integrations for face 3 (x-y plane)
integrations3 = cell(numSubSq3, 1);
for counter = 1 : numSubSq3
    xRange = [minSubX3(counter) maxSubX3(counter)];
    yRange = [minSubY3(counter) maxSubY3(counter)];
    integrations3{counter} = Integrations(xRange, yRange);
end

%% Store results into data structure
squareInfo.midSubX1      = midSubX1;
squareInfo.midSubX2      = midSubX2;
squareInfo.midSubX3      = midSubX3;
squareInfo.midSubY1      = midSubY1;
squareInfo.midSubY2      = midSubY2;
squareInfo.midSubY3      = midSubY3;
squareInfo.midSubSubX1   = midSubSubX1;
squareInfo.midSubSubX2   = midSubSubX2;
squareInfo.midSubSubX3   = midSubSubX3;
squareInfo.midSubSubY1   = midSubSubY1;
squareInfo.midSubSubY2   = midSubSubY2;
squareInfo.midSubSubY3   = midSubSubY3;
squareInfo.minSubSubX1   = minSubSubX1;
squareInfo.minSubSubX2   = minSubSubX2;
squareInfo.minSubSubX3   = minSubSubX3;
squareInfo.maxSubSubX1   = maxSubSubX1;
squareInfo.maxSubSubX2   = maxSubSubX2;
squareInfo.maxSubSubX3   = maxSubSubX3;
squareInfo.minSubSubY1   = minSubSubY1;
squareInfo.minSubSubY2   = minSubSubY2;
squareInfo.minSubSubY3   = minSubSubY3;
squareInfo.maxSubSubY1   = maxSubSubY1;
squareInfo.maxSubSubY2   = maxSubSubY2;
squareInfo.maxSubSubY3   = maxSubSubY3;
squareInfo.surfaceArea1  = surfaceArea1;
squareInfo.surfaceArea2  = surfaceArea2;
squareInfo.surfaceArea3  = surfaceArea3;
squareInfo.integrations1 = integrations1;
squareInfo.integrations2 = integrations2;
squareInfo.integrations3 = integrations3;
squareInfo.numSubSq1     = numSubSq1;
squareInfo.numSubSq2     = numSubSq2;
squareInfo.numSubSq3     = numSubSq3;

end