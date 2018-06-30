% Generate QDFE-center quadrature
function quadrature = QDFEcenter(squareInfo)

% Retrieve needed square properties
minSubSubX   = squareInfo.minSubSubX;
maxSubSubX   = squareInfo.maxSubSubX;
minSubSubY   = squareInfo.minSubSubY;
maxSubSubY   = squareInfo.maxSubSubY;
integrations = squareInfo.integrations;
numSubSq     = squareInfo.numSubSq;

% Initialize storage

xPos    = cell(numSubSq, 1);
yPos    = cell(numSubSq, 1);
gamma   = cell(numSubSq, 3);
theta   = cell(numSubSq, 3);
weights = cell(numSubSq, 1);

% Go through each sub-square
for iSub = 1 : numSubSq
    midX1 = (maxSubSubX{iSub}(1) + minSubSubX{iSub}(1)) / 2;
    midX2 = (maxSubSubX{iSub}(2) + minSubSubX{iSub}(2)) / 2; 
    midX3 = (maxSubSubX{iSub}(3) + minSubSubX{iSub}(3)) / 2; 
    midY1 = (maxSubSubY{iSub}(1) + minSubSubY{iSub}(1)) / 2;
    midY2 = (maxSubSubY{iSub}(4) + minSubSubY{iSub}(4)) / 2; 
    midY3 = (maxSubSubY{iSub}(7) + minSubSubY{iSub}(7)) / 2; 
    xPos{iSub} = [...
        midX1, midX2, midX3, ...
        midX1, midX2, midX3, ...
        midX1, midX2, midX3];
    yPos{iSub} = [...
        midY1, midY1, midY1, ...
        midY2, midY2, midY2, ...
        midY3, midY3, midY3];
    [gamma{iSub, 1}, theta{iSub, 1}] = LocalToGlobal(xPos{iSub}, yPos{iSub});    
    constants = Basis(gamma{iSub, 1}, theta{iSub, 1});     
    weights{iSub, 1} = Weights(integrations{iSub}, constants);         
end

% Rotate to other three faces
for i = 1 : numSubSq
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
quadrature.gamma   = gamma;
quadrature.theta   = theta;
quadrature.weights = weights;

end