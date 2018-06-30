% Second-order spherical harmonics
function intError = SecondInt(quadrature, iRef, intError)

% Retrieve information
gamma = quadrature.gamma;
theta = quadrature.theta;
weights = quadrature.weights;

% Omega-x^2 and omega-x*omega-y integration error
omegaX2 = 0;
omegaXomegaY = 0;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        for k = 1 : 3
            omegaX = cos(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaY = sin(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaX2 = omegaX2 + omegaX ^ 2 * weights{i}(j);
            omegaXomegaY = omegaXomegaY + omegaX * omegaY * weights{i}(j);
        end
    end
end
omegaX2err = abs((pi / 6) - omegaX2) / (pi / 6);
omegaXomegaYerr = abs((1 / 3) - omegaXomegaY) / (1 / 3);

% Store new integration error
if isfield(intError, 'omegaX2err') == 0
    intError.omegaX2err(1) = omegaX2err;
else
    intError.omegaX2err = [intError.omegaX2err, omegaX2err];
end
if isfield(intError, 'omegaXomegaYerr') == 0
    intError.omegaXomegaYerr(1) = omegaXomegaYerr;
else
    intError.omegaXomegaYerr = [intError.omegaXomegaYerr, omegaXomegaYerr];
end

end