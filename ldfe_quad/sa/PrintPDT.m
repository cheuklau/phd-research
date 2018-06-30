% Print quadrature to file for PDT use
function PrintPDT(quadrature, iRef, fid)

% Retrieve information
gamma = quadrature.gamma;
theta = quadrature.theta;
weights = quadrature.weights;

% Normalization factor for PDT
totWgt = 0;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 4
        totWgt = totWgt + weights{i}(j);
    end
end
norm = 0.125 / (3 * totWgt);
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 4
        weights{i}(j) = norm * weights{i}(j);
    end
end

% Quadrature storage for PDT
pdtQuad = zeros((iRef + 1) ^ 2 * 4 * 3, 4);
counter = 1;
for i = 1 : (iRef + 1) ^ 2
    for j = 1 : 4
        for k = 1 : 3
            omegaX = cos(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaY = sin(gamma{i, k}(j)) * sin(theta{i, k}(j));
            omegaZ = cos(theta{i, k}(j));
            pdtQuad(counter, 1) = omegaX;
            pdtQuad(counter, 2) = omegaY;
            pdtQuad(counter, 3) = omegaZ;
            pdtQuad(counter, 4) = weights{i}(j);
            counter = counter + 1;
        end
    end
end

fprintf(fid, '%1.20e %1.20e %1.20e %1.20e\n', pdtQuad');  

end