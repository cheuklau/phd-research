% Generate Quadruple Range quadrature
function quadInfo = QuadRangeQuad()

ref_QR_use = 4;

% QR refinement levels in text files
ref_QR = [2 4 6 8 10 12];

% Number of directions for each refinement level
for i = 1 : size(ref_QR, 2)
    numLines(i) = ref_QR(i) * (ref_QR(i) + 1) / 2;
end

% Import QR data
dataQR = dlmread('QR.DAT');

% Starting line
lineStart = 1;
for j = 1 : ref_QR_use - 1
    lineStart = lineStart + numLines(j);
end

% Store data values
quadInfo = zeros(numLines(ref_QR_use), 5);
weightTotal = 0;
for j = 1 : numLines(ref_QR_use)
    lineUse = lineStart + j - 1;
    quadInfo(j, 1) = 1;
    quadInfo(j, 2) = dataQR(lineUse, 1);
    quadInfo(j, 3) = dataQR(lineUse, 2);
    quadInfo(j, 4) = dataQR(lineUse, 3);
    quadInfo(j, 5) = dataQR(lineUse, 4);
    weightTotal = weightTotal + quadInfo(j, 5);
end

% Normalize weights
for j = 1 : numLines(ref_QR_use)
    quadInfo(j, 5) = quadInfo(j, 5) * (4 * pi / 8) / weightTotal;
end

%{
% Number of angles
nAng = 104;

% Read in quadrature
data = load('QR_26_4.txt');

% Select needed angles
quadInfo = zeros(nAng, 5);
for i = 1 : nAng
    quadInfo(i, 1) = i;
    quadInfo(i, 2) = data(i, 1);
    quadInfo(i, 3) = data(i, 2);
    quadInfo(i, 4) = data(i, 3);
    quadInfo(i, 5) = data(i, 4);
end
%}


end