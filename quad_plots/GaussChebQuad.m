% Generate triangular Gauss-Chebyshev quadrature set
function quadInfo = GaussChebQuad()

ref_GLC_use = 2;

% GLC refinement levels in text files
ref_GLC = [4 8 12 16 32 44 48 64 90 96 128];

% Number of directions for each refinement level
for i = 1 : size(ref_GLC, 2)
    numLines(i) = ref_GLC(i) * (ref_GLC(i) + 1) / 2;
end

% Import GLC data
dataGLC = dlmread('GLC.DAT');

% Starting line
lineStart = 1;
for j = 1 : ref_GLC_use - 1
    lineStart = lineStart + numLines(j);
end

% Store data values
quadInfo = zeros(numLines(ref_GLC_use), 5);
weightTotal = 0;
for j = 1 : numLines(ref_GLC_use)
    lineUse = lineStart + j - 1;
    quadInfo(j, 1) = j;
    quadInfo(j, 2) = dataGLC(lineUse, 1);
    quadInfo(j, 3) = dataGLC(lineUse, 2);
    quadInfo(j, 4) = dataGLC(lineUse, 3);
    quadInfo(j, 5) = dataGLC(lineUse, 4);
    weightTotal = weightTotal + quadInfo(j, 5);
end

% Normalize weights
for j = 1 : numLines(ref_GLC_use)
    quadInfo(j, 5) = quadInfo(j, 5) * (4 * pi / 8) / weightTotal;
end

%{
% Number of polar levels
nPolar = 8;

% Number of angles
nAngles = 0;
for i = 1 : nPolar
    nAngles = nAngles + i;
end

% Resize storage
quadInfo = zeros(nAngles, 5);

% Generate Gauss quadrature for polar angles
[polarPos, polarWgts] = GaussQuad(nPolar, 0, 1);

% Generate azimuthal quadrature
nAzi = nPolar;
dir = 1;
for i = 1 : nPolar
   deltaPhi = pi / (4 * nAzi);
   aziAng = deltaPhi;
   aziWgt = 2*deltaPhi;
   cosTheta = polarPos(i);
   sinTheta = sqrt(1 - cosTheta ^ 2);
   for k = 1 : nAzi
      quadInfo(dir, 1) = dir;
      quadInfo(dir, 2) = cos(aziAng) * sinTheta;
      quadInfo(dir, 3) = sin(aziAng) * sinTheta;
      quadInfo(dir, 4) = cosTheta;
      quadInfo(dir, 5) = aziWgt * polarWgts(i);
      aziAng = aziAng + 2 * deltaPhi;      
      dir = dir + 1;
   end
   nAzi = nAzi - 1;
end
%}

end

