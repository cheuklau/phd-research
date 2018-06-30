% Integrate basis function terms
function integrations = Integrations(xRange, yRange)

% Initialize storage
integrations = zeros(1, 4);

% Gauss quadrature order
N = 32;

% Gauss quadratures
[fx, wx] = lgwt(N, xRange(1), xRange(2));
[fy, wy] = lgwt(N, yRange(1), yRange(2));

% Numerically integrate the basis function terms
counter = 1;
for i = 1 : N
    for j = 1 : N
        temp1(counter) = wx(i) * wy(j) * (1 / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ (3 / 2)));
        temp2(counter) = wx(i) * wy(j) * (1 / (3 * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ 2));
        temp3(counter) = wx(i) * wy(j) * (fx(i) / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ 2));
        temp4(counter) = wx(i) * wy(j) * (fy(j) / (sqrt(3) * (1 / 3 + fx(i) ^ 2 + fy(j) ^ 2) ^ 2));
        counter = counter + 1;
    end
end

% Order temp in increasing order
[~, idx] = sort(abs(temp1));
temp1 = temp1(idx);
[~, idx] = sort(abs(temp2));
temp2 = temp2(idx);
[~, idx] = sort(abs(temp3));
temp3 = temp3(idx);
[~, idx] = sort(abs(temp4));
temp4 = temp4(idx);

for i = 1 : size(temp1, 2)
    integrations(1) = integrations(1) + temp1(i);
    integrations(2) = integrations(2) + temp2(i);
    integrations(3) = integrations(3) + temp3(i);
    integrations(4) = integrations(4) + temp4(i);
end

end