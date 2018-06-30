% Third-order spherical harmonics
function intError = ThirdInt(quadrature, iRef, intError)

% Retrieve information
gamma = quadrature.gamma;
theta = quadrature.theta;
weights = quadrature.weights;

% Omega-x^3, omega-x^2*omega-y and omega-x*omega-y*omega-z integration error
omegaX3 = 0;
omegaX2omegaY = 0;
omegaXomegaYomegaZ = 0;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 9
        for k = 1 : 3
            omegaX = cos(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaY = sin(gamma{i, k}(j)) * sin(theta{i, k}(j)); 
            omegaZ = cos(theta{i, k}(j));
            omegaX3 = omegaX3 + omegaX ^ 3 * weights{i}(j);
            omegaX2omegaY = omegaX2omegaY + omegaX ^ 2 * omegaY * weights{i}(j);
            omegaXomegaYomegaZ = omegaXomegaYomegaZ + omegaX * omegaY * omegaZ * weights{i}(j);
        end        
    end
end
omegaX3err = abs((pi / 8) - omegaX3) / (pi / 8);
omegaX2omegaYerr = abs((pi / 16) - omegaX2omegaY) / (pi / 16);
omegaXomegaYomegaZerr = abs((1 / 8) - omegaXomegaYomegaZ) / (1 / 8);

% Store new integration error
if isfield(intError, 'omegaX3err') == 0
    intError.omegaX3err(1) = omegaX3err;
else
    intError.omegaX3err = [intError.omegaX3err, omegaX3err];
end
if isfield(intError, 'omegaX2omegaYerr') == 0
    intError.omegaX2omegaYerr(1) = omegaX2omegaYerr;
else
    intError.omegaX2omegaYerr = [intError.omegaX2omegaYerr, omegaX2omegaYerr];
end
if isfield(intError, 'omegaXomegaYomegaZerr') == 0
    intError.omegaXomegaYomegaZerr(1) = omegaXomegaYomegaZerr;
else
    intError.omegaXomegaYomegaZerr = [intError.omegaXomegaYomegaZerr, omegaXomegaYomegaZerr];
end

end