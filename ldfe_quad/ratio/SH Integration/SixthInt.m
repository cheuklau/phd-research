% Sixth-order spherical harmonics
function intError = SixthInt(quadrature, iRef, intError)

% Retrieve information
gamma = quadrature.gamma;
theta = quadrature.theta;
weights = quadrature.weights;

% Omega-x^6, omega-x^5*omega-y, omega-x^4*omega-y^2,
% omega-x^4*omega-y*omega-z, omega-x^3*omega-y^2*omega-z,
% omega-x^3*omega-y^2*omega-z and omega-x^2*omega-y^2*omega-z^2
% integration error
omegaX6 = 0;
omegaX5omegaY = 0;
omegaX4omegaY2 = 0;
omegaX4omegaYomegaZ = 0;
omegaX3omegaY2omegaZ = 0;
omegaX2omegaY2omegaZ2 = 0;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 4
        for k = 1 : 3
            omegaX = cos(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaY = sin(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaZ = cos(theta{i, k}(j));
            omegaX6 = omegaX6 + omegaX ^ 6 * weights{i}(j);
            omegaX5omegaY = omegaX5omegaY + omegaX ^ 5 * omegaY * weights{i}(j);
            omegaX4omegaY2 = omegaX4omegaY2 + omegaX ^ 4 * omegaY ^ 2 * weights{i}(j);
            omegaX4omegaYomegaZ = omegaX4omegaYomegaZ + omegaX ^ 4 * omegaY * omegaZ * weights{i}(j);
            omegaX3omegaY2omegaZ = omegaX3omegaY2omegaZ + omegaX ^ 3 * omegaY ^ 2 * omegaZ * weights{i}(j);
            omegaX2omegaY2omegaZ2 = omegaX2omegaY2omegaZ2 + omegaX ^ 2 * omegaY ^ 2 * omegaZ ^ 2 * weights{i}(j);
        end        
    end
end
omegaX6err = abs((pi / 14) - omegaX6) / (pi / 14);
omegaX5omegaYerr = abs((8 / 105) - omegaX5omegaY) / (8 / 105);
omegaX4omegaY2err = abs((pi / 70) - omegaX4omegaY2) / (pi / 70);
omegaX4omegaYomegaZerr = abs((1 / 35) - omegaX4omegaYomegaZ) / (1 / 35);
omegaX3omegaY2omegaZerr = abs((2 / 105) - omegaX3omegaY2omegaZ) / (2 / 105);
omegaX2omegaY2omegaZ2err = abs((pi / 210) - omegaX2omegaY2omegaZ2) / (pi / 210);

% Store new integration error
if isfield(intError, 'omegaX6err') == 0
    intError.omegaX6err(1) = omegaX6err;
else
    intError.omegaX6err = [intError.omegaX6err, omegaX6err];
end
if isfield(intError, 'omegaX5omegaYerr') == 0
    intError.omegaX5omegaYerr(1) = omegaX5omegaYerr;
else
    intError.omegaX5omegaYerr = [intError.omegaX5omegaYerr, omegaX5omegaYerr];
end
if isfield(intError, 'omegaX4omegaY2err') == 0
    intError.omegaX4omegaY2err(1) = omegaX4omegaY2err;
else
    intError.omegaX4omegaY2err = [intError.omegaX4omegaY2err, omegaX4omegaY2err];
end
if isfield(intError, 'omegaX4omegaYomegaZerr') == 0
    intError.omegaX4omegaYomegaZerr(1) = omegaX4omegaYomegaZerr;
else
    intError.omegaX4omegaYomegaZerr = [intError.omegaX4omegaYomegaZerr, omegaX4omegaYomegaZerr];
end
if isfield(intError, 'omegaX3omegaY2omegaZerr') == 0
    intError.omegaX3omegaY2omegaZerr(1) = omegaX3omegaY2omegaZerr;
else
    intError.omegaX3omegaY2omegaZerr = [intError.omegaX3omegaY2omegaZerr, omegaX3omegaY2omegaZerr];
end
if isfield(intError, 'omegaX2omegaY2omegaZ2err') == 0
    intError.omegaX2omegaY2omegaZ2err(1) = omegaX2omegaY2omegaZ2err;
else
    intError.omegaX2omegaY2omegaZ2err = [intError.omegaX2omegaY2omegaZ2err, omegaX2omegaY2omegaZ2err];
end

end