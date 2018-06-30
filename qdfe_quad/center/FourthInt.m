% Fourth-order spherical harmonics
function intError = FourthInt(quadrature, iRef, intError)

% Retrieve information
gamma = quadrature.gamma;
theta = quadrature.theta;
weights = quadrature.weights;

% Omega-x^4, omega-x^3*omega-y, omega-x^2*omega-y^2 and
% omega-x^2*omega-y*omega-z integration error
omegaX4 = 0; 
omegaX3omegaY = 0;
omegaX2omegaY2 = 0;
omegaX2omegaYomegaZ = 0;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        for k = 1 : 3
            omegaX = cos(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaY = sin(gamma{i, k}(j)) * sin(theta{i, k}(j)); 
            omegaZ = cos(theta{i, k}(j));
            omegaX4 = omegaX4 + omegaX ^ 4 * weights{i}(j);
            omegaX3omegaY = omegaX3omegaY + omegaX ^ 3 * omegaY * weights{i}(j);
            omegaX2omegaY2 = omegaX2omegaY2 + omegaX ^ 2 * omegaY ^ 2 * weights{i}(j);
            omegaX2omegaYomegaZ = omegaX2omegaYomegaZ + omegaX ^ 2 * omegaY * omegaZ * weights{i}(j);
        end
    end
end
omegaX4err = abs((pi / 10) - omegaX4) / (pi / 10);
omegaX3omegaYerr = abs((2 / 15) - omegaX3omegaY) / (2 / 15);
omegaX2omegaY2err = abs((pi / 30) - omegaX2omegaY2) / (pi / 30);
omegaX2omegaYomegaZerr = abs((1 / 15) - omegaX2omegaYomegaZ) / (1 / 15);

% Store new integration error
if isfield(intError, 'omegaX4err') == 0
    intError.omegaX4err(1) = omegaX4err;
else
    intError.omegaX4err = [intError.omegaX4err, omegaX4err];
end
if isfield(intError, 'omegaX3omegaYerr') == 0
    intError.omegaX3omegaYerr(1) = omegaX3omegaYerr;
else
    intError.omegaX3omegaYerr = [intError.omegaX3omegaYerr, omegaX3omegaYerr];
end
if isfield(intError, 'omegaX2omegaY2err') == 0
    intError.omegaX2omegaY2err(1) = omegaX2omegaY2err;
else
    intError.omegaX2omegaY2err = [intError.omegaX2omegaY2err, omegaX2omegaY2err];
end
if isfield(intError, 'omegaX2omegaYomegaZerr') == 0
    intError.omegaX2omegaYomegaZerr(1) = omegaX2omegaYomegaZerr;
else
    intError.omegaX2omegaYomegaZerr = [intError.omegaX2omegaYomegaZerr, omegaX2omegaYomegaZerr];
end

end