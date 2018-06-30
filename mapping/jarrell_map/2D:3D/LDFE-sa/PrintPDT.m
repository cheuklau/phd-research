% Print quadrature to file for PDT use
function PrintPDT(quadrature, order, fid)

% Retrieve information
gamma    = quadrature.gamma;
theta    = quadrature.theta;
weights  = quadrature.weights;
basis    = quadrature.basis;
numSubSq = (order + 1) ^ 2;

% Normalization factor for PDT
totWgt = 0;
for i = 1 : numSubSq
    for j = 1 : 4
        totWgt = totWgt + weights{i}(j);
    end
end
norm = 0.125 / (3 * totWgt);
for i = 1 : numSubSq
    for j = 1 : 4
        weights{i}(j) = norm * weights{i}(j);
    end
end

% Quadrature storage
pdtQuad = zeros(numSubSq * 4 * 3, 8);
counter = 1;
for iFace = 1 : 3
    for iSub = 1 : numSubSq
        for iDir = 1 : 4
            omegaX = cos(gamma{iSub, iFace}(iDir)) * sin(theta{iSub, iFace}(iDir));
            omegaY = sin(gamma{iSub, iFace}(iDir)) * sin(theta{iSub, iFace}(iDir));
            omegaZ = cos(theta{iSub, iFace}(iDir));
            pdtQuad(counter, 1) = omegaX;
            pdtQuad(counter, 2) = omegaY;
            pdtQuad(counter, 3) = omegaZ;
            pdtQuad(counter, 4) = weights{iSub}(iDir);
            % Basis function constants must be rotated since faces are
            % rotated
            if iFace == 1
                pdtQuad(counter, 5) = basis{iSub}(1, iDir);
                pdtQuad(counter, 6) = basis{iSub}(2, iDir);
                pdtQuad(counter, 7) = basis{iSub}(3, iDir);
                pdtQuad(counter, 8) = basis{iSub}(4, iDir);
            elseif iFace == 2
                pdtQuad(counter, 5) = basis{iSub}(1, iDir);
                pdtQuad(counter, 6) = basis{iSub}(4, iDir);
                pdtQuad(counter, 7) = basis{iSub}(2, iDir);
                pdtQuad(counter, 8) = basis{iSub}(3, iDir);
            elseif iFace == 3
                pdtQuad(counter, 5) = basis{iSub}(1, iDir);
                pdtQuad(counter, 6) = basis{iSub}(3, iDir);
                pdtQuad(counter, 7) = basis{iSub}(4, iDir);
                pdtQuad(counter, 8) = basis{iSub}(2, iDir);
            end
                
            counter = counter + 1;
        end
    end
end

fprintf(fid, '%1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e \n', pdtQuad');  

end