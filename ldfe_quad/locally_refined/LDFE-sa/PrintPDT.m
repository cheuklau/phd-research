% Print quadrature to file for PDT use
function PrintPDT(squareInfo, quadrature, fid)

%% Retrieve quadrature data
gamma1    = quadrature.gamma1;
gamma2    = quadrature.gamma2;
gamma3    = quadrature.gamma3;
theta1    = quadrature.theta1;
theta2    = quadrature.theta2;
theta3    = quadrature.theta3;
weights1  = quadrature.weights1;
weights2  = quadrature.weights2;
weights3  = quadrature.weights3;
numSubSq1 = squareInfo.numSubSq1;
numSubSq2 = squareInfo.numSubSq2;
numSubSq3 = squareInfo.numSubSq3;

%% Normalization factor for PDT
totWgt = 0;
for i = 1 : numSubSq1
    for j = 1 : 4
        totWgt = totWgt + weights1{i}(j);
    end
end
for i = 1 : numSubSq2
    for j = 1 : 4
        totWgt = totWgt + weights2{i}(j);
    end
end
for i = 1 : numSubSq3
    for j = 1 : 4
        totWgt = totWgt + weights3{i}(j);
    end
end
norm = 1 / totWgt;

%% Normalize weights
for i = 1 : numSubSq1
    for j = 1 : 4
        weights1{i}(j) = norm * weights1{i}(j);
    end
end
for i = 1 : numSubSq2
    for j = 1 : 4
        weights2{i}(j) = norm * weights2{i}(j);
    end
end
for i = 1 : numSubSq3
    for j = 1 : 4
        weights3{i}(j) = norm * weights3{i}(j);
    end
end

%% Quadrature storage for PDT
pdtQuad = zeros(numSubSq1 + numSubSq2 + numSubSq3, 4);
counter = 1;
for i = 1 : numSubSq1
    for j = 1 : 4
        omegaX = cos(gamma1{i}(j)) * sin(theta1{i}(j));
        omegaY = sin(gamma1{i}(j)) * sin(theta1{i}(j));
        omegaZ = cos(theta1{i}(j));
        pdtQuad(counter, 1) = omegaX;
        pdtQuad(counter, 2) = omegaY;
        pdtQuad(counter, 3) = omegaZ;
        pdtQuad(counter, 4) = weights1{i}(j);
        counter = counter + 1;
    end
end
for i = 1 : numSubSq2
    for j = 1 : 4
        omegaX = cos(gamma2{i}(j)) * sin(theta2{i}(j));
        omegaY = sin(gamma2{i}(j)) * sin(theta2{i}(j));
        omegaZ = cos(theta2{i}(j));
        pdtQuad(counter, 1) = omegaX;
        pdtQuad(counter, 2) = omegaY;
        pdtQuad(counter, 3) = omegaZ;
        pdtQuad(counter, 4) = weights2{i}(j);
        counter = counter + 1;
    end
end
for i = 1 : numSubSq3
    for j = 1 : 4
        omegaX = cos(gamma3{i}(j)) * sin(theta3{i}(j));
        omegaY = sin(gamma3{i}(j)) * sin(theta3{i}(j));
        omegaZ = cos(theta3{i}(j));
        pdtQuad(counter, 1) = omegaX;
        pdtQuad(counter, 2) = omegaY;
        pdtQuad(counter, 3) = omegaZ;
        pdtQuad(counter, 4) = weights3{i}(j);
        counter = counter + 1;
    end
end

% Display the total number of directions
fprintf('the total number of sub-squares per octant: %i \n', numSubSq1 + numSubSq2 + numSubSq3);

fprintf(fid, '%1.20e %1.20e %1.20e %1.20e\n', pdtQuad');

end