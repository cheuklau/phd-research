% First-order spherical harmonics
function intError = FirstInt(quadrature, iRef, intError)

% Retrieve information
gamma = quadrature.gamma;
theta = quadrature.theta;
weights = quadrature.weights;

% Omega-x integration error
omegaX = 0;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        for k = 1 : 3
            omegaX = omegaX + cos(gamma{i, k}(j)) * sin(theta{i, k}(j)) * weights{i}(j);
        end
    end
end
omegaXerr = abs((pi / 4) - omegaX) / (pi / 4);

% Store new integration error
if isfield(intError, 'omegaXerr') == 0
    intError.omegaXerr(1) = omegaXerr;
else
    intError.omegaXerr = [intError.omegaXerr, omegaXerr];
end

end