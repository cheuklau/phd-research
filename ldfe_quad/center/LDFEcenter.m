% Generate LDFE-center quadrature
function quadrature = LDFEcenter(squareInfo)

% Retrieve needed square properties
midSubSubX   = squareInfo.midSubSubX;
midSubSubY   = squareInfo.midSubSubY;
integrations = squareInfo.integrations;
numSubSq     = squareInfo.numSubSq;

% Initialize storage
xPos    = cell(numSubSq, 1);
yPos    = cell(numSubSq, 1);
gamma   = cell(numSubSq, 3);
theta   = cell(numSubSq, 3);
weights = cell(numSubSq, 1);

% Go through each sub-square
for i = 1 : numSubSq
    xPos{i} = [midSubSubX{i}(1), midSubSubX{i}(2), midSubSubX{i}(3), midSubSubX{i}(4)]; 
    yPos{i} = [midSubSubY{i}(1), midSubSubY{i}(2), midSubSubY{i}(3), midSubSubY{i}(4)]; 
    [gamma{i, 1}, theta{i, 1}] = LocalToGlobal(xPos{i}, yPos{i});
    constants = Basis(gamma{i, 1}, theta{i, 1});    
    weights{i, 1} = Weights(integrations{i}, constants);
end

% Rotate to other three faces
for i = 1 : numSubSq
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