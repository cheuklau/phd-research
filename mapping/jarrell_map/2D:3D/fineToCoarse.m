function solCoarse = fineToCoarse(solFrom, orderFrom, orderTo, mapping)

% Fine quadrature information
orderFine = orderFrom;
quadFine  = GenQuad(orderFine);
solFine   = solFrom;

% Coarse quadrature information
orderCoarse = orderFine - 1;
quadCoarse = GenQuad(orderCoarse);

% Map coarser one level at a time
while orderCoarse >= orderTo
    
    % Fine quadrature information for a single face
    numSubSqFine = size(quadFine, 1) / (3 * 4);
    numSubSqFineRow = numSubSqFine ^ 0.5;
    
    % Coarse quadrature information for a single face
    numSubSqCoarse = size(quadCoarse, 1) / (3 * 4);
    
    % Initialize storage for mapped solution
    solCoarse = cell(3, 1);
    
    % Go through each face
    for iFace = 1 : 3
        
        % Go through each coarse sub-square of current face
        deltaX  = 0;
        deltaY  = 0;
        for i = 1 : numSubSqCoarse
            
            % Gather quadrature data for current coarse sub-square
            omegaXcoarse = zeros(4, 1);
            omegaYcoarse = zeros(4, 1);
            omegaZcoarse = zeros(4, 1);
            weightcoarse = zeros(4, 1);
            for j = 1 : 4
                iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                omegaXcoarse(j) = quadCoarse{iDir}(1);
                omegaYcoarse(j) = quadCoarse{iDir}(2);
                omegaZcoarse(j) = quadCoarse{iDir}(3);
                weightcoarse(j) = quadCoarse{iDir}(4);
            end
            
            % Fine quadrature data lines corresponding to coarse quadrature
            lineFine(1) = (iFace - 1) * numSubSqFine * 4 + deltaX * 4 + deltaY * numSubSqFineRow * 4 + 1;
            lineFine(2) = lineFine(1) + 4;
            lineFine(3) = lineFine(1) + numSubSqFineRow * 4;
            lineFine(4) = lineFine(3) + 4;
            
            % Preserve 0th and 1st moments
            if mapping == 1
                
                % A-matrix
                aMat = zeros(4);
                for j = 1 : 4
                    aMat(1, j) = weightcoarse(j);
                    aMat(2, j) = weightcoarse(j) * omegaXcoarse(j);
                    aMat(3, j) = weightcoarse(j) * omegaYcoarse(j);
                    aMat(4, j) = weightcoarse(j) * omegaZcoarse(j);
                end
                
                % b-vector
                bvec = zeros(4, 1);
                for j = 1 : 4
                    for k = 0 : 3
                        lineUse = lineFine(j) + k;
                        bvec(1) = bvec(1) + quadFine{lineUse}(4) * solFine(lineUse);
                        bvec(2) = bvec(2) + quadFine{lineUse}(4) * quadFine{lineUse}(1) * solFine(lineUse);
                        bvec(3) = bvec(3) + quadFine{lineUse}(4) * quadFine{lineUse}(2) * solFine(lineUse);
                        bvec(4) = bvec(4) + quadFine{lineUse}(4) * quadFine{lineUse}(3) * solFine(lineUse);
                    end
                end
                
                % Solve for mapped coarse values
                solCoarse{iFace}{i} = aMat \ bvec;
                
                % Preserve 1st moment using LDFE basis functions
            elseif mapping == 2
                
                % Interpolate fine LDFE basis functions for coarse solution
                solCoarse{iFace}{i} = zeros(4, 1);
                for j = 1 : 4
                    for k = 0 : 3
                        lineUse = lineFine(j) + k;
                        solCoarse{iFace}{i}(j) = ...
                            solCoarse{iFace}{i}(j) + ...
                            solFine(lineUse) * ( ...
                            quadFine{lineUse}(5) + ...
                            quadFine{lineUse}(6) * omegaXcoarse(j) + ...
                            quadFine{lineUse}(7) * omegaYcoarse(j) + ...
                            quadFine{lineUse}(8) * omegaZcoarse(j));
                    end
                end
                
                % Conservation factor numerator
                numerator = 0;
                for j = 1 : 4
                    for k = 0 : 3
                        lineUse = lineFine(j) + k;
                        numerator = ...
                            numerator + ...
                            quadFine{lineUse}(1) * ...
                            solFine(lineUse) * quadFine{lineUse}(4);
                    end
                end
                
                % Conservation factor denominator
                denominator = 0;
                for j = 1 : 4
                    denominator = ...
                        denominator + ...
                        omegaXcoarse(j) * ...
                        solCoarse{iFace}{i}(j) * weightcoarse(j);
                end
                
                % Apply conservation factor to mapped solution               
                if abs(denominator) ~= 0
                    factor = numerator / denominator;
                    solCoarse{iFace}{i} = solCoarse{iFace}{i} * factor;
                end
                
            end
            
            % Update x and y positions
            deltaX = deltaX + 2;
            if deltaX == numSubSqFineRow
                deltaX = 0;
                deltaY = deltaY + 2;
            end
            
        end
        
    end
    
    % Update coarse and fine information
    if orderCoarse > orderTo
        
        % Fine quadrature order
        orderFine = orderFine - 1;
        
        % Fine quadrature
        quadFine = GenQuad(orderFine);
        
        % Fine quadrature solution
        % Note: Have to rearrange to match solFrom format
        numSubSqFine = size(quadFine, 1) / (3 * 4);
        counter = 1;
        for i = 1 : 3
            for j = 1 : numSubSqFine
                for k = 1 : 4
                    solFine(counter) = solCoarse{i}{j}(k);
                    counter = counter + 1;
                end
            end
        end
        
        % Coarse quadrature order
        orderCoarse = orderCoarse - 1;
        
        % Coarse quadrature
        quadCoarse = GenQuad(orderCoarse);
        
    else
        
        % Coarse quadrature order
        orderCoarse = orderCoarse - 1;
        
    end
    
end

end