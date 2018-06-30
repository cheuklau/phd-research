% Plot selected quadrature over the first octant
% Reads in PDT generated quadrature in the format:
% <angle number> <omega-x> <omega-y> <omega-z> <weight>

% Read directional cosines and weight from file
%quadInfo = dlmread('LSquad.txt');
%quadInfo = GaussChebQuad();
%quadInfo = QuadRangeQuad();
%quadInfo = QuadRangeQuad_BAILEY();
quadInfo = OctLDFE();
%quadInfo = OctLDFEcube();
%quadInfo = OctQDFEcube();
%quadInfo = OctLDFEcubeLocal();
%quadInfo = GLLDFE();
%quadInfo = OctPGLC();
numAngles = size(quadInfo, 1);
xPos = zeros(1, numAngles);
yPos = zeros(1, numAngles);
zPos = zeros(1, numAngles);
weight = zeros(1, numAngles);
for i = 1 : 1 : numAngles
   % Calculate theta
   theta = acos(quadInfo(i, 4));
   % Calculate gamma
   gamma = asin(quadInfo(i, 3) / sin(theta));
   % Convert from spherical to Cartesian coordinates for plotting
   [xPos(i), yPos(i), zPos(i)] = sph2cart(gamma, (pi / 2 - theta), 1); 
   % Store the weights
   weight(i) = quadInfo(i, 5);
end

% Plot the results
figure;
scatter3(xPos, yPos, zPos, weight * 90000, 'k', 'Filled');
hold on;
set(gca, 'FontSize', 16);
xlabel('\mu','FontSize',24);
ylabel('\eta','FontSize',24);
zlabel('\xi','FontSize',24);