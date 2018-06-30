% Print quadrature to file for PDT use
function PrintPDT(quadrature, order, fid)

% Retrieve information
gamma    = quadrature.gamma;
theta    = quadrature.theta;
weights  = quadrature.weights;
basis    = quadrature.basis;
basis2   = quadrature.basis2;
basis3   = quadrature.basis3;
numSubSq = (order + 1) ^ 2;

% Normalization factor for PDT
totWgt = 0;
for i = 1 : numSubSq
    for j = 1 : 9
        totWgt = totWgt + weights{i}(j);
    end
end
norm = 0.125 / (3 * totWgt);
for i = 1 : numSubSq
    for j = 1 : 9
        weights{i}(j) = norm * weights{i}(j);
    end
end

% Quadrature storage
pdtQuad = zeros(numSubSq * 9 * 3, 13);
counter = 1;
for iFace = 1: 3
    for iSub = 1 : numSubSq
        for iDir = 1 : 9
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
                pdtQuad(counter, 5)  = basis{iSub}(1, iDir); % 1
                pdtQuad(counter, 6)  = basis{iSub}(2, iDir); % x
                pdtQuad(counter, 7)  = basis{iSub}(3, iDir); % y
                pdtQuad(counter, 8)  = basis{iSub}(4, iDir); % z
                pdtQuad(counter, 9)  = basis{iSub}(5, iDir); % x^2-y^2
                pdtQuad(counter, 10) = basis{iSub}(6, iDir); % z^2
                pdtQuad(counter, 11) = basis{iSub}(7, iDir); % x*y
                pdtQuad(counter, 12) = basis{iSub}(8, iDir); % x*z
                pdtQuad(counter, 13) = basis{iSub}(9, iDir); % y*z                
            elseif iFace == 2
                % Rotations:
                % x becomes z
                % y becomes x
                % z becomes y
                pdtQuad(counter, 5)  = basis2{iSub}(1, iDir); % 1
                pdtQuad(counter, 6)  = basis2{iSub}(2, iDir); % z
                pdtQuad(counter, 7)  = basis2{iSub}(3, iDir); % x
                pdtQuad(counter, 8)  = basis2{iSub}(4, iDir); % y
                pdtQuad(counter, 9)  = basis2{iSub}(5, iDir); % z^2-x^2
                pdtQuad(counter, 10) = basis2{iSub}(6, iDir); % y^2
                pdtQuad(counter, 11) = basis2{iSub}(7, iDir); % x*z
                pdtQuad(counter, 12) = basis2{iSub}(8, iDir); % y*z
                pdtQuad(counter, 13) = basis2{iSub}(9, iDir); % x*y
            elseif iFace == 3
                % Rotations:
                % x becomes y
                % y becomes z
                % z becomes x
                pdtQuad(counter, 5)  = basis3{iSub}(1, iDir); % 1
                pdtQuad(counter, 6)  = basis3{iSub}(2, iDir); % y
                pdtQuad(counter, 7)  = basis3{iSub}(3, iDir); % z
                pdtQuad(counter, 8)  = basis3{iSub}(4, iDir); % x
                pdtQuad(counter, 9)  = basis3{iSub}(5, iDir); % y^2-z^2
                pdtQuad(counter, 10) = basis3{iSub}(6, iDir); % x^2
                pdtQuad(counter, 11) = basis3{iSub}(7, iDir); % y*z
                pdtQuad(counter, 12) = basis3{iSub}(8, iDir); % x*y
                pdtQuad(counter, 13) = basis3{iSub}(9, iDir); % x*z 
            end
            
            counter = counter + 1;
            
        end
    end
end

fprintf(fid, '%1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e %1.20e\n', pdtQuad');  

end