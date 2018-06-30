% Basis functions for a sub-square
function [basis, basis2, basis3] = Basis(gamma, theta)

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

% Test -- basis functions for second face
% x becomes z
% y becomes x
% z becomes y
A = zeros(9, 9);
for i = 1 : 9
    x = cos(gamma(i)) * sin(theta(i));
    y = sin(gamma(i)) * sin(theta(i));
    z = cos(theta(i));
    A(i, 1) = 1;
    A(i, 2) = z;
    A(i, 3) = x;
    A(i, 4) = y;
    A(i, 5) = z ^ 2 - x ^ 2;
    A(i, 6) = y ^ 2;
    A(i, 7) = x * z;
    A(i, 8) = y * z;
    A(i, 9) = x * y;
end

b = eye(9);

if rcond(A) < 1e-12 || isnan(rcond(A)) == 1
    basis2 = zeros(9);
else
    basis2 = A \ b;
end

% Test -- basis functions for third face    
% x becomes y
% y becomes z
% z becomes x
A = zeros(9, 9);
for i = 1 : 9
    x = cos(gamma(i)) * sin(theta(i));
    y = sin(gamma(i)) * sin(theta(i));
    z = cos(theta(i));
    A(i, 1) = 1;
    A(i, 2) = y;
    A(i, 3) = z;
    A(i, 4) = x;
    A(i, 5) = y ^ 2 - z ^ 2;
    A(i, 6) = x ^ 2;
    A(i, 7) = y * z;
    A(i, 8) = x * y;
    A(i, 9) = x * z;
end

b = eye(9);

if rcond(A) < 1e-12 || isnan(rcond(A)) == 1
    basis3 = zeros(9);
else
    basis3 = A \ b;
end

end