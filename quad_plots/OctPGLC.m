% Product GLC
function quadInfo = OctPGLC()

% Number of angles
nAng = 24;

% Read in quadrature
data = load('PGLCREC.txt');

% Select needed angles
quadInfo = zeros(nAng, 5);
tot_weight = 0;
for i = 1 : nAng
    quadInfo(i, 1) = i;
    quadInfo(i, 2) = data(i, 1);
    quadInfo(i, 3) = data(i, 2);
    quadInfo(i, 4) = data(i, 3);
    quadInfo(i, 5) = data(i, 4);
    tot_weight = tot_weight + quadInfo(i, 5);
end

% Normalize weights
factor = 1 / tot_weight;
for i = 1 : nAng
    quadInfo(i, 5) = data(i, 4) * factor;
end


end