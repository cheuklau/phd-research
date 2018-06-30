% Calculate residual
function res = calc_res(...
    ratio_corner, ratio_side, ...
    subSubLengthX, subSubLengthY, ...
    minSubSubX, maxSubSubX, ...
    minSubSubY, maxSubSubY, ...
    integrations, surfaceArea)

% x and y positions on sub-square
xPos = [...
    maxSubSubX(1) - subSubLengthX(1) * ratio_corner, ...
   (maxSubSubX(2) + minSubSubX(2)) / 2, ...
    minSubSubX(3) + subSubLengthX(3) * ratio_corner, ...
    ...
    maxSubSubX(4) - subSubLengthX(4) * ratio_side, ...
   (maxSubSubX(5) + minSubSubX(5)) / 2, ...
    minSubSubX(6) + subSubLengthX(6) * ratio_side, ...
    ...
    maxSubSubX(7) - subSubLengthX(7) * ratio_corner, ...
   (maxSubSubX(8) + minSubSubX(8)) / 2, ...
    minSubSubX(9) + subSubLengthX(9) * ratio_corner];

yPos = [...
    maxSubSubY(1) - subSubLengthY(1) * ratio_corner, ...
    maxSubSubY(2) - subSubLengthY(2) * ratio_side, ...
    maxSubSubY(3) - subSubLengthY(3) * ratio_corner, ...
    ...
   (maxSubSubY(4) + minSubSubY(4)) / 2, ...
   (maxSubSubY(5) + minSubSubY(5)) / 2, ...
   (maxSubSubY(6) + minSubSubY(6)) / 2, ...
    ...
    minSubSubY(7) + subSubLengthY(7) * ratio_corner, ...
    minSubSubY(8) + subSubLengthY(8) * ratio_side, ...
    minSubSubY(9) + subSubLengthY(9) * ratio_corner];

% Global angles
[gamma, theta] = LocalToGlobal(xPos, yPos);

% Basis functions
constants = Basis(gamma, theta);

% Weights
weights = Weights(integrations, constants);

% Residual of center sub-sub-square
res = (surfaceArea(5) - weights(5)) / surfaceArea(5);

end