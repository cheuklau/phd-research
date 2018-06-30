% Basis functions for a sub-square
function basis = Basis(gamma, theta)

% Create matrix
A = zeros(4, 4);
for i = 1 : 4
    A(i, 1) = 1;
    A(i, 2) = cos(gamma(i)) * sin(theta(i));
    A(i, 3) = sin(gamma(i)) * sin(theta(i));
    A(i, 4) = cos(theta(i));
end

b = eye(4);

if rcond(A) < 1e-12 || isnan(rcond(A)) == 1
    basis = zeros(4);
else
    basis = A \ b;
end

end