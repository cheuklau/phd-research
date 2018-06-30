function squareInfo = SquareProp(iRef)

%% Basic properties

% Square length
length = 1 / sqrt(3);

% Total number of sub-squares
numSubSq = (iRef + 1) ^ 2;

% Number of sub-sub-squares
numSubSubSq = 9 * numSubSq;

% Number of rings
numRings = sqrt(numSubSubSq);

% Ideal weight of each ordinate
wgtIdeal = pi / (6 * numSubSubSq);

% Iteration criteria
delta = 1e-12;
maxIter = 100;

%% Generate ring limits
counterSub     = 2; % Keeps track of sub-square limits
counterSubSub  = 1; % Keeps track of sub-square middle positions
subSqLimits    = zeros(iRef + 2, 1);       % Stores sub-square limits
subSqLimits(1) = length;                   % Outermost sub-square position
subSubSqLimits = zeros(2 * (iRef + 1), 1); % Stores sub-sub-square limits

% Go through each ring
for i = 1 : numRings - 1
    
    % Reset ratio limits, convergence tracker and iteration counter   
    upperLimit = 1.0;
    lowerLimit = 0;
    converged  = 0;
    counter    = 0;
    
    % Ideal weight for current plus previous rings
    wgtIdealTotal = wgtIdeal * i ^ 2;
    
    % Start bisection method
    while converged == 0 && counter < maxIter
        
        % Calculate new ratio
        ratio = (upperLimit + lowerLimit) / 2; 
        
        % Calculate SQ surface area
        surfaceArea = calcSA(ratio * length, length, ratio * length, length);
        
        % Residual        
        residual = (surfaceArea - wgtIdealTotal) / wgtIdealTotal;
                
        % Check for convergence
        if abs(residual) < delta
            
           converged = 1; 
            
           % Store as sub-sub-square limit if ring number not div by 3
           if mod(i, 3) ~= 0
               
               subSubSqLimits(counterSubSub) = ratio * length;
               
               counterSubSub = counterSubSub + 1;
               
               % Store as sub-square limit for even number ring
           else
               
               subSqLimits(counterSub) = ratio * length;
               
               counterSub = counterSub + 1;
               
           end
           
        end
        
        % Surface area too large (increase ratio)
        if (residual > 0)
            
            lowerLimit = ratio;
            
                                              
            % Surface area too small (decrease ratio)
        else
            
            upperLimit = ratio;
                                    
        end
        
        % Increment iteration counter
        counter = counter + 1;                
        
    end        

end

% Reverse the order of the sub-square middle position and limits
subSubSqLimits = flipud(subSubSqLimits);
subSqLimits    = flipud(subSqLimits);

% Store sub-square limits
minSubX = zeros(numSubSq, 1);
maxSubX = zeros(numSubSq, 1);
minSubY = zeros(numSubSq, 1);
maxSubY = zeros(numSubSq, 1);
counter = 1;
for i = 1 : iRef + 1
    for j = 1 : iRef + 1
        minSubX(counter) = subSqLimits(j);
        maxSubX(counter) = subSqLimits(j + 1);
        minSubY(counter) = subSqLimits(i);
        maxSubY(counter) = subSqLimits(i + 1);
        counter = counter + 1;
    end
end

% Sub-sub-square limits
minSubSubX = cell(numSubSq, 1);
maxSubSubX = cell(numSubSq, 1);
minSubSubY = cell(numSubSq, 1);
maxSubSubY = cell(numSubSq, 1);
subSubLengthX = cell(numSubSq, 1);
subSubLengthY = cell(numSubSq, 1);
counter = 1;
for i = 1 : iRef + 1    
    for j = 1 : iRef + 1
        x1 = subSqLimits(j);
        x2 = subSubSqLimits((j - 1) * 2 + 1);
        x3 = subSubSqLimits((j - 1) * 2 + 2);
        x4 = subSqLimits(j + 1);
        y1 = subSqLimits(i);
        y2 = subSubSqLimits((i - 1) * 2 + 1);
        y3 = subSubSqLimits((i - 1) * 2 + 2);
        y4 = subSqLimits(i + 1);
        % Lower-left sub-sub-square
        minSubSubX{counter}(1) = x1;
        maxSubSubX{counter}(1) = x2;
        minSubSubY{counter}(1) = y1;
        maxSubSubY{counter}(1) = y2;
        subSubLengthX{counter}(1) = x2 - x1;
        subSubLengthY{counter}(1) = y2 - y1;
        % Lower-center sub-sub-square
        minSubSubX{counter}(2) = x2;
        maxSubSubX{counter}(2) = x3;
        minSubSubY{counter}(2) = y1;
        maxSubSubY{counter}(2) = y2;
        subSubLengthX{counter}(2) = x3 - x2;
        subSubLengthY{counter}(2) = y2 - y1;
        % Lower-right sub-sub-square
        minSubSubX{counter}(3) = x3;
        maxSubSubX{counter}(3) = x4;
        minSubSubY{counter}(3) = y1;
        maxSubSubY{counter}(3) = y2;
        subSubLengthX{counter}(3) = x4 - x3;
        subSubLengthY{counter}(3) = y2 - y1;
        % Center-left sub-sub-square
        minSubSubX{counter}(4) = x1;
        maxSubSubX{counter}(4) = x2;
        minSubSubY{counter}(4) = y2;
        maxSubSubY{counter}(4) = y3;
        subSubLengthX{counter}(4) = x2 - x1;
        subSubLengthY{counter}(4) = y3 - y2;
        % Center-center sub-sub-square
        minSubSubX{counter}(5) = x2;
        maxSubSubX{counter}(5) = x3;    
        minSubSubY{counter}(5) = y2;    
        maxSubSubY{counter}(5) = y3;
        subSubLengthX{counter}(5) = x3 - x2;
        subSubLengthY{counter}(5) = y3 - y2;
        % Center-right sub-sub-square
        minSubSubX{counter}(6) = x3;   
        maxSubSubX{counter}(6) = x4;    
        minSubSubY{counter}(6) = y2;    
        maxSubSubY{counter}(6) = y3;
        subSubLengthX{counter}(6) = x4 - x3;
        subSubLengthY{counter}(6) = y3 - y2;
        % Upper-left sub-sub-square
        minSubSubX{counter}(7) = x1;    
        maxSubSubX{counter}(7) = x2;    
        minSubSubY{counter}(7) = y3;    
        maxSubSubY{counter}(7) = y4;
        subSubLengthX{counter}(7) = x2 - x1;
        subSubLengthY{counter}(7) = y4 - y3;
        % Upper-center sub-sub-square
        minSubSubX{counter}(8) = x2;    
        maxSubSubX{counter}(8) = x3;    
        minSubSubY{counter}(8) = y3;    
        maxSubSubY{counter}(8) = y4; 
        subSubLengthX{counter}(8) = x3 - x2;
        subSubLengthY{counter}(8) = y4 - y3;
        % Upper-right sub-sub-square
        minSubSubX{counter}(9) = x3;    
        maxSubSubX{counter}(9) = x4;
        minSubSubY{counter}(9) = y3;
        maxSubSubY{counter}(9) = y4;
        subSubLengthX{counter}(9) = x4 - x3;
        subSubLengthY{counter}(9) = y4 - y3;
        counter = counter + 1;
    end
end
          
% Sub-sub-square surface area
N = 16;
SA_max = 0;
SA_min = 999;
surfaceArea = cell(numSubSq, 1);
for k = 1 : numSubSq
    surfaceArea{k} = zeros(9, 1);
    for m = 1 : 9
        [fx, wx] = lgwt(N, minSubSubX{k}(m), maxSubSubX{k}(m));
        [fy, wy] = lgwt(N, minSubSubY{k}(m), maxSubSubY{k}(m));
        for i = 1 : N
            for j = 1 : N
                temp = (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2)))...
                    * wx(i) * wy(j);
                surfaceArea{k}(m) = surfaceArea{k}(m) + temp;
            end
        end
        if surfaceArea{k}(m) > SA_max
            SA_max = surfaceArea{k}(m);
        elseif surfaceArea{k}(m) < SA_min
            SA_min = surfaceArea{k}(m);
        end
    end
end

fprintf('the max-to-min ratio for sub-sub-square areas: %f \n', SA_max/SA_min);

% Basis function integration over sub-square
integrations = cell(numSubSq, 1);
for counter = 1 : numSubSq
    xRange = [minSubX(counter) maxSubX(counter)];
    yRange = [minSubY(counter) maxSubY(counter)];
    integrations{counter} = Integrations(xRange, yRange);
end

% Store results into data structure
squareInfo.minSubSubX    = minSubSubX;
squareInfo.maxSubSubX    = maxSubSubX;
squareInfo.minSubSubY    = minSubSubY;
squareInfo.maxSubSubY    = maxSubSubY;
squareInfo.subSubLengthX = subSubLengthX;
squareInfo.subSubLengthY = subSubLengthY;
squareInfo.surfaceArea   = surfaceArea;
squareInfo.integrations  = integrations;
squareInfo.numSubSq      = numSubSq;

end