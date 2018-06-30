% Generate 2D RZ GL-LDFE quadrature
function quadInfo = GLLDFE()

% Polar level selection
polarIndex = 2;

% Number of polar levels
numPolar = [2, 4, 8, 16, 32, 64];

% Number of angles
nPolar = numPolar(polarIndex);
nAngles = nPolar * (nPolar + 1);

% Starting line of data
lineStart = 1;
for i = 1 : polarIndex - 1
    lineStart = lineStart + numPolar(i) * (numPolar(i) + 1);    
end

% Read in quadrature
data = load('GLLDFE.dat');

% Select needed angles
quadInfo = zeros(nAngles, 5);
for i = 1 : nAngles
    row = lineStart + i - 1;
    quadInfo(i, 1) = i;
    quadInfo(i, 2) = data(row, 1);
    quadInfo(i, 3) = data(row, 2);
    quadInfo(i, 4) = data(row, 3);
    quadInfo(i, 5) = data(row, 4);
end

end

