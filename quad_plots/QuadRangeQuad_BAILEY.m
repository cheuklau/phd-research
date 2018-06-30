% Generate Quadruple Range quadrature
function quadInfo = QuadRangeQuad_BAILEY()

% Number of angles
nAng = 72;

% Read in quadrature
data = load('QR_18_4.txt');

% Select needed angles
quadInfo = zeros(nAng, 5);
for i = 1 : nAng
    quadInfo(i, 1) = i;
    quadInfo(i, 2) = data(i, 1);
    quadInfo(i, 3) = data(i, 2);
    quadInfo(i, 4) = data(i, 3);
    quadInfo(i, 5) = data(i, 4);
end



end