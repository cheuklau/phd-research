function [minSubX, maxSubX, minSubY, maxSubY, ...
          midSubX, midSubY, midSubSubX, midSubSubY, ...
          minSubSubX, maxSubSubX, minSubSubY, maxSubSubY, ...
          xPos, yPos, ...
          numSubSq] = calcUniform(iRef)
    
% Square length
length = 1 / sqrt(3);
 
% Number of sub-squares
numSubSq = (iRef + 1) ^ 2;

% Number of sub-sub-squares
numSubSubSq = 4 * numSubSq;

% Number of rings
numRings = sqrt(numSubSubSq);

% Ideal weight of each ordinate
wgtIdeal = pi / (6 * numSubSubSq);

% Iteration criteria
delta   = 1e-12;
maxIter = 100;

% Generate ring limits
counterMidSub  = 1; % Keeps track of sub-square middle positions
counterSub     = 2; % Keeps track of sub-square limits
subSqLimits    = zeros(iRef + 2, 1); % Stores sub-square limits
subSqLimits(1) = length;             % Outermost position
midSubSq       = zeros(iRef + 1, 1); % Stores sub-square middle positions

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
            
            % Store as sub-square middle position for odd number ring
            if mod(i, 2) ~= 0
                
                midSubSq(counterMidSub) = ratio * length;
                
                counterMidSub = counterMidSub + 1;
                
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
midSubSq    = flipud(midSubSq);
subSqLimits = flipud(subSqLimits);

% Store sub-square limits
minSubX = zeros(numSubSq, 1);
maxSubX = zeros(numSubSq, 1);
minSubY = zeros(numSubSq, 1);
maxSubY = zeros(numSubSq, 1);
midSubX = zeros(numSubSq, 1);
midSubY = zeros(numSubSq, 1);
xPos   = zeros(numSubSq, 1);
yPos   = zeros(numSubSq, 1);
counter = 1;
for i = 1 : iRef + 1
    for j = 1 : iRef + 1
        minSubX(counter) = subSqLimits(j);
        maxSubX(counter) = subSqLimits(j + 1);
        midSubX(counter) = midSubSq(j);
        minSubY(counter) = subSqLimits(i);
        maxSubY(counter) = subSqLimits(i + 1);
        midSubY(counter) = midSubSq(i);
        xPos(counter) = j;
        yPos(counter) = i;
        counter = counter + 1;
    end
end

% Store sub-sub-square midpoints and limits
midSubSubX = cell(numSubSq, 1);
midSubSubY = cell(numSubSq, 1);
minSubSubX = cell(numSubSq, 1);
maxSubSubX = cell(numSubSq, 1);
minSubSubY = cell(numSubSq, 1);
maxSubSubY = cell(numSubSq, 1);
for i = 1 : numSubSq
    
    minSubSubX{i} = [midSubX(i), midSubX(i), minSubX(i), minSubX(i)];
    maxSubSubX{i} = [maxSubX(i), maxSubX(i), midSubX(i), midSubX(i)];
    minSubSubY{i} = [midSubY(i), minSubY(i), minSubY(i), midSubY(i)];
    maxSubSubY{i} = [maxSubY(i), midSubY(i), midSubY(i), maxSubY(i)];
    
    midSubSubX{i} = [...
        (maxSubSubX{i}(1) + minSubSubX{i}(1)) / 2, ...
        (maxSubSubX{i}(2) + minSubSubX{i}(2)) / 2, ...
        (maxSubSubX{i}(3) + minSubSubX{i}(3)) / 2, ...
        (maxSubSubX{i}(4) + minSubSubX{i}(4)) / 2];
    midSubSubY{i} = [...
        (maxSubSubY{i}(1) + minSubSubY{i}(1)) / 2, ...
        (maxSubSubY{i}(2) + minSubSubY{i}(2)) / 2, ...
        (maxSubSubY{i}(3) + minSubSubY{i}(3)) / 2, ...
        (maxSubSubY{i}(4) + minSubSubY{i}(4)) / 2];
    
end

end