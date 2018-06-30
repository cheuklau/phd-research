% Fifth-order spherical harmonics
function intError = FifthInt(quadrature, iRef, intError)

% Retrieve information
gamma = quadrature.gamma;
theta = quadrature.theta;
weights = quadrature.weights;

% Omega-x^5, omega-x^4*omega-y, omega-x^3*omega-y^2,
% omega-x^3*omega-y*omega-z and omega-x^2*omega-y^2*omega-z^2
% integration error
omegaX5 = 0;
omegaX4omegaY = 0;
omegaX3omegaY2 = 0;
omegaX3omegaYomegaZ = 0;
omegaX2omegaY2omegaZ = 0;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        for k = 1 : 3
            omegaX = cos(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaY = sin(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaZ = cos(theta{i, k}(j));
            omegaX5 = omegaX5 + omegaX ^ 5 * weights{i}(j);
            omegaX4omegaY = omegaX4omegaY + omegaX ^ 4 * omegaY * weights{i}(j);
            omegaX3omegaY2 = omegaX3omegaY2 + omegaX ^ 3 * omegaY ^ 2 * weights{i}(j);
            omegaX3omegaYomegaZ = omegaX3omegaYomegaZ + omegaX ^ 3 * omegaY * omegaZ * weights{i}(j);
            omegaX2omegaY2omegaZ = omegaX2omegaY2omegaZ + omegaX ^ 2 * omegaY ^ 2 * omegaZ * weights{i}(j);
        end
    end
end
omegaX5err = abs((pi / 12) - omegaX5) / (pi / 12);
omegaX4omegaYerr = abs((pi / 32) - omegaX4omegaY) / (pi / 32);
omegaX3omegaY2err = abs((pi / 48) - omegaX3omegaY2) / (pi / 48);
omegaX3omegaYomegaZerr = abs((1 / 24) - omegaX3omegaYomegaZ) / (1 / 24);
omegaX2omegaY2omegaZerr = abs((pi / 96) - omegaX2omegaY2omegaZ) / (pi / 96);

% Store new integration error
if isfield(intError, 'omegaX5err') == 0
    intError.omegaX5err(1) = omegaX5err;
else
    intError.omegaX5err = [intError.omegaX5err, omegaX5err];
end
if isfield(intError, 'omegaX4omegaYerr') == 0
    intError.omegaX4omegaYerr(1) = omegaX4omegaYerr;
else
    intError.omegaX4omegaYerr = [intError.omegaX4omegaYerr, omegaX4omegaYerr];
end
if isfield(intError, 'omegaX3omegaY2err') == 0
    intError.omegaX3omegaY2err(1) = omegaX3omegaY2err;
else
    intError.omegaX3omegaY2err = [intError.omegaX3omegaY2err, omegaX3omegaY2err];
end
if isfield(intError, 'omegaX3omegaYomegaZerr') == 0
    intError.omegaX3omegaYomegaZerr(1) = omegaX3omegaYomegaZerr;
else
    intError.omegaX3omegaYomegaZerr = [intError.omegaX3omegaYomegaZerr, omegaX3omegaYomegaZerr];
end
if isfield(intError, 'omegaX2omegaY2omegaZerr') == 0
    intError.omegaX2omegaY2omegaZerr(1) = omegaX2omegaY2omegaZerr;
else
    intError.omegaX2omegaY2omegaZerr = [intError.omegaX2omegaY2omegaZerr, omegaX2omegaY2omegaZerr];
end

end