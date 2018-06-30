% Generate octahedron-based LDFE quadratures
function quadInfo = OctLDFEcube()

% Refinement level
refLevel = 5;

% Number of angles
nAng = (refLevel + 1) ^ 2 * 4 * 3;

% Starting line of data
lineStart = 1;
for i = 0 : 1 : refLevel - 1
    lineStart = lineStart + (i + 1) ^ 2 * 4 * 3;
end

% Read in quadrature
data = load('LDFEsa0to5L.DAT');

% Select needed angles
quadInfo = zeros(nAng, 5);
for i = 1 : nAng
    row = lineStart + i - 1;
    quadInfo(i, 1) = i;
    quadInfo(i, 2) = data(row, 1);
    quadInfo(i, 3) = data(row, 2);
    quadInfo(i, 4) = data(row, 3);
    quadInfo(i, 5) = data(row, 4);
end

end
