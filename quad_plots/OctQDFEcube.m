% Generate octahedron-based QDFE quadratures
function quadInfo = OctQDFEcube()

% Refinement level
refLevel = 3;

% Number of angles
nAng = (refLevel + 1) ^ 2 * 9 * 3;

% Starting line of data
lineStart = 1;
for i = 0 : refLevel - 1
    lineStart = lineStart + (i + 1) ^ 2 * 9 * 3;
end

% Read in quadrature
data = load('QDFEsaEven.DAT');

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
