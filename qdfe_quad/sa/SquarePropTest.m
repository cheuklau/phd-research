% Calculate basic properties of the square
function squareInfo = SquarePropTest(iRef)

% Square length
length = 1 / sqrt(3);

%{
% Sub-square length
subLength = length / (iRef + 1);

% Even sub-square limits
limits = zeros(iRef+2,1);
for i = 2 : iRef + 2
    limits(i) = limits(i - 1) + subLength;
end
counter = 1;
minSubX = zeros((iRef + 1) ^ 2, 1);
maxSubX = zeros((iRef + 1) ^ 2, 1);
minSubY = zeros((iRef + 1) ^ 2, 1);
maxSubY = zeros((iRef + 1) ^ 2, 1);
for i = 1 : iRef + 1
    for j = 1 : iRef + 1
        minSubX(counter) = limits(j);
        maxSubX(counter) = limits(j + 1);
        minSubY(counter) = limits(i);
        maxSubY(counter) = limits(i + 1);
        counter = counter + 1;
    end
end

% Sub-square limits forcing diagonal SQ surface areas to be equal
[minSubX, maxSubX, minSubY, maxSubY] = SquarePropEven(iRef, minSubX, maxSubX, minSubY, maxSubY);
%}

% Even arc-length sub-square limits
delta_angle = (pi / 4) / (iRef + 1);
current_angle = 0;
limits = zeros(iRef + 2, 1);
for i = 1 : iRef + 2
    limits(i) = length * tan(current_angle);
    current_angle = current_angle + delta_angle;
end
counter = 1;
minSubX = zeros((iRef + 1) ^ 2, 1);
maxSubX = zeros((iRef + 1) ^ 2, 1);
minSubY = zeros((iRef + 1) ^ 2, 1);
maxSubY = zeros((iRef + 1) ^ 2, 1);
for i = 1 : iRef + 1
    for j = 1 : iRef + 1
        minSubX(counter) = limits(j);
        maxSubX(counter) = limits(j + 1);
        minSubY(counter) = limits(i);
        maxSubY(counter) = limits(i + 1);
        counter = counter + 1;
    end
end

% Sub-sub-square limits
minSubSubX = cell((iRef + 1) ^ 2, 1);
maxSubSubX = cell((iRef + 1) ^ 2, 1);
minSubSubY = cell((iRef + 1) ^ 2, 1);
maxSubSubY = cell((iRef + 1) ^ 2, 1);
for counter = 1 : (iRef + 1) ^ 2
    subSubLengthX = (maxSubX(counter) - minSubX(counter)) / 3;
    subSubLengthY = (maxSubY(counter) - minSubY(counter)) / 3;
    minSubSubX{counter}(1) = minSubX(counter);    
    maxSubSubX{counter}(1) = minSubX(counter) + subSubLengthX;    
    minSubSubY{counter}(1) = minSubY(counter);    
    maxSubSubY{counter}(1) = minSubY(counter) + subSubLengthY;  
    minSubSubX{counter}(2) = minSubX(counter) + subSubLengthX;    
    maxSubSubX{counter}(2) = minSubX(counter) + 2 * subSubLengthX;    
    minSubSubY{counter}(2) = minSubY(counter);    
    maxSubSubY{counter}(2) = minSubY(counter) + subSubLengthY;
    minSubSubX{counter}(3) = minSubX(counter) + 2 * subSubLengthX;    
    maxSubSubX{counter}(3) = minSubX(counter) + 3 * subSubLengthX;    
    minSubSubY{counter}(3) = minSubY(counter);    
    maxSubSubY{counter}(3) = minSubY(counter) + subSubLengthY;    
    minSubSubX{counter}(4) = minSubX(counter);
    maxSubSubX{counter}(4) = minSubX(counter) + subSubLengthX;
    minSubSubY{counter}(4) = minSubY(counter) + subSubLengthY;
    maxSubSubY{counter}(4) = minSubY(counter) + 2 * subSubLengthY;
    minSubSubX{counter}(5) = minSubX(counter) + subSubLengthX;    
    maxSubSubX{counter}(5) = minSubX(counter) + 2 * subSubLengthX;    
    minSubSubY{counter}(5) = minSubY(counter) + subSubLengthY;    
    maxSubSubY{counter}(5) = minSubY(counter) + 2 * subSubLengthY;    
    minSubSubX{counter}(6) = minSubX(counter) + 2 * subSubLengthX;   
    maxSubSubX{counter}(6) = minSubX(counter) + 3 * subSubLengthX;    
    minSubSubY{counter}(6) = minSubY(counter) + subSubLengthY;    
    maxSubSubY{counter}(6) = minSubY(counter) + 2 * subSubLengthY;    
    minSubSubX{counter}(7) = minSubX(counter);    
    maxSubSubX{counter}(7) = minSubX(counter) + subSubLengthX;    
    minSubSubY{counter}(7) = minSubY(counter) + 2 * subSubLengthY;    
    maxSubSubY{counter}(7) = minSubY(counter) + 3 * subSubLengthY;    
    minSubSubX{counter}(8) = minSubX(counter) + subSubLengthX;    
    maxSubSubX{counter}(8) = minSubX(counter) + 2 * subSubLengthX;    
    minSubSubY{counter}(8) = minSubY(counter) + 2 * subSubLengthY;    
    maxSubSubY{counter}(8) = minSubY(counter) + 3 * subSubLengthY;    
    minSubSubX{counter}(9) = minSubX(counter) + 2 * subSubLengthX;    
    maxSubSubX{counter}(9) = minSubX(counter) + 3 * subSubLengthX;
    minSubSubY{counter}(9) = minSubY(counter) + 2 * subSubLengthY;
    maxSubSubY{counter}(9) = minSubY(counter) + 3 * subSubLengthY;
end

% Sub-sub-square midpoints
midSubSubX = cell((iRef + 1) ^ 2, 1);
midSubSubY = cell((iRef + 1) ^ 2, 1);
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        midSubSubX{i}(j) = (maxSubSubX{i}(j) + minSubSubX{i}(j)) / 2;
        midSubSubY{i}(j) = (maxSubSubY{i}(j) + minSubSubY{i}(j)) / 2;
    end
end

% Sub-sub-square surface area
N = 16;
SA_max = 0;
SA_min = 999;
surfaceArea = cell((iRef + 1) ^ 2, 1);
for k = 1 : (iRef + 1) ^ 2
    surfaceArea{k} = zeros(9, 1);
    SA_sub = 0;
    for m = 1 : 9
        [fx, wx] = lgwt(N, minSubSubX{k}(m), maxSubSubX{k}(m));
        [fy, wy] = lgwt(N, minSubSubY{k}(m), maxSubSubY{k}(m));
        for i = 1 : N
            for j = 1 : N
                temp = (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2))) ...
                    * wx(i) * wy(j);
                surfaceArea{k}(m) = surfaceArea{k}(m) + temp;
            end
        end
        SA_sub = SA_sub + surfaceArea{k}(m);
    end
    if SA_sub > SA_max
        SA_max = SA_sub;
    elseif SA_sub < SA_min
        SA_min = SA_sub;
    end    
end

fprintf('the max-to-min ratio for sub-square areas: %f \n', SA_max/SA_min);

% Basis function integration over sub-square
integrations = cell((iRef + 1) ^ 2, 1);
for counter = 1 : (iRef + 1) ^ 2
    xRange = [minSubX(counter) maxSubX(counter)];
    yRange = [minSubY(counter) maxSubY(counter)];
    integrations{counter} = Integrations(xRange, yRange);
end

% Store results into data structure
squareInfo.length = length;
squareInfo.minSubX = minSubX;
squareInfo.maxSubX = maxSubX;
squareInfo.minSubY = minSubY;
squareInfo.maxSubY = maxSubY;
squareInfo.minSubSubX = minSubSubX;
squareInfo.maxSubSubX = maxSubSubX;
squareInfo.minSubSubY = minSubSubY;
squareInfo.maxSubSubY = maxSubSubY;
squareInfo.midSubSubX = midSubSubX;
squareInfo.midSubSubY = midSubSubY;
squareInfo.surfaceArea = surfaceArea;
squareInfo.integrations = integrations;

end