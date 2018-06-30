function quad = GenQuad(order)

% Number of LDFE ranges
numRanges = 2 ^ (order + 1);

% Angles for each LDFE range
deltaMu = (2 / numRanges) / 4;
muCurrent = -1;
for i = 1 : numRanges
    mu{i}(1) = muCurrent + deltaMu;
    mu{i}(2) = mu{i}(1) + 2 * deltaMu;
    muCurrent = muCurrent + 4 * deltaMu;
end

% Basis functions for each LDFE range
aMat = zeros(2);
bMat = eye(2);
for i = 1 : numRanges    
    
    % Assemble A-matrix
    for j = 1 : 2
        aMat(j, 1) = 1;        
        aMat(j, 2) = mu{i}(j);
    end
    
    % Solve for basis function constants
    basis{i} = aMat \ bMat;        
    
    % Mu range
    maxMu = mu{i}(2) + deltaMu;
    minMu = mu{i}(1) - deltaMu;
    
    % Solve for the weights
    for j = 1 : 2        
        weight{i}(j) = basis{i}(1, j) * (maxMu - minMu) + ...
            (basis{i}(2, j) / 2) * (maxMu ^ 2 - minMu ^ 2);        
    end
    
end

% Store solutions
quad.mu = mu;
quad.basis = basis;
quad.weight = weight;

end