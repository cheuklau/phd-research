% Integrate basis function terms
function integrations = Integrations(xRange, yRange)

% Initialize storage
integrations = zeros(1, 4);

% Gauss quadrature order
N = 64;

% Gauss quadratures
[fx, wx] = lgwt(N, xRange(1), xRange(2));
[fy, wy] = lgwt(N, yRange(1), yRange(2));

% Numerically integrate the basis function terms
for i = 1 : N
    for j = 1 : N
        integrations(1) = integrations(1) + wx(i) * wy(j) * (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2)));
        integrations(2) = integrations(2) + wx(i) * wy(j) * (1 / (3 * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ 2));
        integrations(3) = integrations(3) + wx(i) * wy(j) * (fx(i) / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ 2));
        integrations(4) = integrations(4) + wx(i) * wy(j) * (fy(j) / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ 2));
    end
end

end