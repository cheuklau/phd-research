% Generate octahedron-based LDFE quadratures
function quadInfo = OctLDFE()

% Refinement level
refLevel = 2;

% Number of angles
nAng = 4 ^ (refLevel + 1);

% Starting line of data
lineStart = 1;
for i = 0 : refLevel - 1
    lineStart = lineStart + 4 ^ (i + 1);
end

% Read in quadrature
data = load('LDFE.dat');

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
