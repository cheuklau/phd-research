% Basis functions for a sub-square
function basis = Basis(gamma, theta)

% Create matrix
A = zeros(9, 9);
for i = 1 : 9
    x = cos(gamma(i)) * sin(theta(i));
    y = sin(gamma(i)) * sin(theta(i));
    z = cos(theta(i));
    A(i, 1) = 1;
    A(i, 2) = x;
    A(i, 3) = y;
    A(i, 4) = z;
    A(i, 5) = x ^ 2 - y ^ 2;
    A(i, 6) = z ^ 2;
    A(i, 7) = x * y;
    A(i, 8) = x * z;
    A(i, 9) = y * z;
end

b = eye(9);

if rcond(A) < 1e-12 || isnan(rcond(A)) == 1
    basis = zeros(9);
else
    basis = A \ b;
end

end