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

basis = A \ b;

end