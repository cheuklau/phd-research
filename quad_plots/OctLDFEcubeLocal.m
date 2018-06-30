% Generate octahedron-based LDFE quadratures
function quadInfo = OctLDFEcubeLocal()

% Refinement to plot
ref = 4;

% Number of sub-squares
numSubSq = [24 36 78 165]; % +60% feathered

numDirs = numSubSq * 4;

% Starting line of data
lineStart = 1;
for i = 1 : ref - 1
    lineStart = lineStart + numDirs(i);
end

% Read in quadrature
data = load('LDFEQUADRATIOLOCAL60FEATHER.DAT');

% Select needed angles
quadInfo = zeros(numDirs(ref), 5);
for i = 1 : numDirs(ref)
    row = lineStart + i - 1;
    quadInfo(i, 1) = i;
    quadInfo(i, 2) = data(row, 1);
    quadInfo(i, 3) = data(row, 2);
    quadInfo(i, 4) = data(row, 3);
    quadInfo(i, 5) = data(row, 4);
end

end
