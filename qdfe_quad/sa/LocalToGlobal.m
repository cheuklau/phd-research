% Convert local (x,y) to global angles
function [gamma, theta] = LocalToGlobal(xPos, yPos)

% Initialize storage
gamma = zeros(1, 9);
theta = zeros(1, 9);

% Convert local (x,y) to global angles
for i = 1 : 9
    gamma(i) = atan(sqrt(3) * xPos(i));     
    theta(i) = (pi / 2) - atan(sqrt(3) * yPos(i) / sqrt(1 + 3 * xPos(i) ^ 2));
end

end