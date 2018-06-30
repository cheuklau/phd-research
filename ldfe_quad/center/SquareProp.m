% Calculate basic properties of the square
function squareInfo = SquareProp(iRef)

% Square length
length = 1 / sqrt(3);

% Sub-square length
subLength = length / (iRef + 1);

% Total number of sub-squares
numSubSq = (iRef + 1) ^ 2;

% Sub-square limits
limits = zeros(iRef + 2, 1);
for i = 2 : iRef + 2
    limits(i) = limits(i - 1) + subLength;
end
counter = 1;
minSubX = zeros(numSubSq, 1);
maxSubX = zeros(numSubSq, 1);
minSubY = zeros(numSubSq, 1);
maxSubY = zeros(numSubSq, 1);
for i = 1 : iRef + 1
    for j = 1 : iRef + 1
        minSubX(counter) = limits(j);
        maxSubX(counter) = limits(j + 1);
        minSubY(counter) = limits(i);
        maxSubY(counter) = limits(i + 1);
        counter = counter + 1;
    end
end

% Sub-square midpoints
mid = zeros((iRef + 1) * 2 + 1);
for i = 2 : (iRef + 1) * 2 + 1
    mid(i) = mid(i - 1) + subLength / 2;
end
counter = 1;
midSubX = zeros(numSubSq, 1);
midSubY = zeros(numSubSq, 1);
for i = 1 : iRef + 1
    for j = 1 : iRef + 1
        midSubX(counter) = mid(2 * j);
        midSubY(counter) = mid(2 * i);
        counter = counter + 1;
    end
end

% Sub-sub-square midpoints
midSubSubX = cell(numSubSq, 1);
midSubSubY = cell(numSubSq, 1);
for counter = 1 : numSubSq
    midSubSubX{counter} = [midSubX(counter) + subLength / 4, ...
        midSubX(counter) + subLength / 4, ...
        midSubX(counter) - subLength / 4, ...
        midSubX(counter) - subLength / 4];
    midSubSubY{counter} = [midSubY(counter) + subLength / 4, ...
        midSubY(counter) - subLength / 4, ...
        midSubY(counter) - subLength / 4, ...
        midSubY(counter) + subLength / 4];
end

% Sub-sub-square limits
minSubSubX = cell(numSubSq, 1);
maxSubSubX = cell(numSubSq, 1);
minSubSubY = cell(numSubSq, 1);
maxSubSubY = cell(numSubSq, 1);
for counter = 1 : numSubSq
    minSubSubX{counter} = [midSubSubX{counter}(1) - subLength / 4, ...
        midSubSubX{counter}(2) - subLength / 4, ...
        midSubSubX{counter}(3) - subLength / 4, ...
        midSubSubX{counter}(4) - subLength / 4];
    maxSubSubX{counter} = [midSubSubX{counter}(1) + subLength / 4, ...
        midSubSubX{counter}(2) + subLength / 4, ...
        midSubSubX{counter}(3) + subLength / 4, ...
        midSubSubX{counter}(4) + subLength / 4];
    minSubSubY{counter} = [midSubSubY{counter}(1) - subLength / 4,...
        midSubSubY{counter}(2) - subLength / 4,...
        midSubSubY{counter}(3) - subLength / 4,...
        midSubSubY{counter}(4) - subLength / 4];
    maxSubSubY{counter} = [midSubSubY{counter}(1) + subLength / 4, ...
        midSubSubY{counter}(2) + subLength / 4, ...
        midSubSubY{counter}(3) + subLength / 4, ...
        midSubSubY{counter}(4) + subLength / 4];
end

% Sub-sub-square surface area
N = 64;
surfaceArea = cell(numSubSq, 1);
for k = 1 : numSubSq
    surfaceArea{k} = zeros(4, 1);
    for m = 1 : 4
        [fx, wx] = lgwt(N, minSubSubX{k}(m), maxSubSubX{k}(m));
        [fy, wy] = lgwt(N, minSubSubY{k}(m), maxSubSubY{k}(m));
        for i = 1 : N
            for j = 1 : N
                temp = (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2)))...
                    * wx(i) * wy(j);
                surfaceArea{k}(m) = surfaceArea{k}(m) + temp;
            end
        end
    end
end

% Basis function integration over sub-square
integrations = cell(numSubSq, 1);
for counter = 1 : numSubSq
    xRange = [minSubX(counter) maxSubX(counter)];
    yRange = [minSubY(counter) maxSubY(counter)];
    integrations{counter} = Integrations(xRange, yRange);
end

% Store results into data structure
squareInfo.length = length;
squareInfo.subLength = subLength;
squareInfo.minSubX = minSubX;
squareInfo.maxSubX = maxSubX;
squareInfo.minSubY = minSubY;
squareInfo.maxSubY = maxSubY;
squareInfo.midSubX = midSubX;
squareInfo.midSubY = midSubY;
squareInfo.midSubSubX = midSubSubX;
squareInfo.midSubSubY = midSubSubY;
squareInfo.surfaceArea = surfaceArea;
squareInfo.integrations = integrations;
squareInfo.numSubSq = numSubSq;

end