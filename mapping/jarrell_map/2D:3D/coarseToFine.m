function solFine = coarseToFine(solFrom, orderFrom, orderTo, mapping)

% Coarse quadrature information
orderCoarse = orderFrom;
quadCoarse  = GenQuad(orderCoarse);
solCoarse   = solFrom;

% Fine quadrature information
orderFine = orderFrom + 1;
quadFine  = GenQuad(orderFine);

% While fine order less than quadrature order you are mapping to
while orderFine <= orderTo
    
    % Fine quadrature information for a single face
    numSubSqFine    = size(quadFine, 1) / (3 * 4);
    numSubSqFineRow = numSubSqFine ^ 0.5;
    
    % Coarse quadrature information for a single face
    numSubSqCoarse = size(quadCoarse, 1) / (3 * 4);
    
    % Initialize storage for mapped solution
    solFine = cell(3, 1);
    
    % Go through each face
    for iFace = 1 : 3
        
        % Go through each coarse sub-square of current face
        deltaX = 0;
        deltaY = 0;
        for i = 1 : numSubSqCoarse
            
            % Fine sub-squares corresponding to coarse sub-square
            subSqFine(3) = deltaX + deltaY * numSubSqFineRow + 1;
            subSqFine(2) = subSqFine(3) + 1;
            subSqFine(4) = subSqFine(3) + numSubSqFineRow;
            subSqFine(1) = subSqFine(4) + 1;
            
            % Fine quadrature lines corresponding to coarse quadrature
            lineFine(3) = 4 * ((iFace - 1) * numSubSqFine + deltaX + deltaY * numSubSqFineRow) + 1;
            lineFine(2) = lineFine(3) + 4;
            lineFine(4) = lineFine(3) + numSubSqFineRow * 4;
            lineFine(1) = lineFine(4) + 4;
            
            % Preserve 0th and 1st moments
            if mapping == 1
                
                % Go through each corresponding fine sub-square
                for j = 1 : 4
                    
                    % Corresponding coarse quadrature direction information
                    iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                    omegaXcoarse = quadCoarse{iDir}(1);
                    omegaYcoarse = quadCoarse{iDir}(2);
                    omegaZcoarse = quadCoarse{iDir}(3);
                    weightcoarse = quadCoarse{iDir}(4);
                    solCoarseUse = solCoarse(iDir);
                    
                    % A-matrix
                    aMat = zeros(4);
                    for k = 1 : 4
                        lineUse = lineFine(j) + k - 1;
                        aMat(1, k) = quadFine{lineUse}(4);
                        aMat(2, k) = quadFine{lineUse}(4) * quadFine{lineUse}(1);
                        aMat(3, k) = quadFine{lineUse}(4) * quadFine{lineUse}(2);
                        aMat(4, k) = quadFine{lineUse}(4) * quadFine{lineUse}(3);
                    end
                    
                    % b-vector
                    bvec = zeros(4, 1);
                    bvec(1) = weightcoarse * solCoarseUse;
                    bvec(2) = weightcoarse * omegaXcoarse * solCoarseUse;
                    bvec(3) = weightcoarse * omegaYcoarse * solCoarseUse;
                    bvec(4) = weightcoarse * omegaZcoarse * solCoarseUse;
                    
                    % Solve for mapped coarse values
                    solFine{iFace}{subSqFine(j)} = aMat \ bvec;
                    
                end
                
                % Preserve 1st moment using LDFE basis functions
            elseif mapping == 2
                
                % Initialize conservation factor denominator
                denominator = 0;
                
                % Go through each corresponding fine sub-square
                for j = 1 : 4
                    
                    % Initialize fine solution
                    solFine{iFace}{subSqFine(j)} = zeros(4, 1);
                    
                    % Go through each fine quadrature direction
                    for k = 1 : 4
                        
                        % Go through each coarse quadrature direction
                        lineUse = lineFine(j) + k - 1;
                        for m = 1 : 4
                            iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + m;
                            solFine{iFace}{subSqFine(j)}(k) = ...
                                solFine{iFace}{subSqFine(j)}(k) + ...
                                solCoarse(iDir) * ( ...
                                quadCoarse{iDir}(5) + ...
                                quadCoarse{iDir}(6) * quadFine{lineUse}(1) + ...
                                quadCoarse{iDir}(7) * quadFine{lineUse}(2) + ...
                                quadCoarse{iDir}(8) * quadFine{lineUse}(3));
                        end
                        
                        % Increment denominator
                        denominator = denominator + ...
                            solFine{iFace}{subSqFine(j)}(k) * ...
                            quadFine{lineUse}(4) * ...
                            quadFine{lineUse}(1);
                        
                    end
                    
                end
                
                % Conservation factor numerator
                numerator = 0;
                for j = 1 : 4
                    iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                    numerator = ...
                        numerator + ...
                        solCoarse(iDir) * ...
                        quadCoarse{iDir}(4) * ...
                        quadCoarse{iDir}(1);
                end
                                
                % Apply conservation factor to mapped solution
                if abs(denominator) ~= 0
                    factor = numerator / denominator;
                    for j = 1 : 4
                        solFine{iFace}{subSqFine(j)} = solFine{iFace}{subSqFine(j)} * factor;                        
                    end                                       
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
    if orderFine < orderTo
        
        % Coarse quadrature order
        orderCoarse = orderCoarse + 1;
        
        % Coarse quadrature
        quadCoarse = GenQuad(orderCoarse);
        
        % Coarse quadrature solution
        % Note: Have to rearrange to match solFrom format
        numSubSqCoarse = size(quadCoarse, 1) / (3 * 4);
        counter = 1;
        for i = 1 : 3
            for j = 1 : numSubSqCoarse
                for k = 1 : 4
                    solCoarse(counter) = solFine{i}{j}(k);
                    counter = counter + 1;
                end
            end
        end
        
        % Fine quadrature order
        orderFine = orderFine + 1;
        
        % Fine quadrature
        quadFine = GenQuad(orderFine);
        
    else
        
        % Fine quadrature order
        orderFine = orderFine + 1;
        
    end
    
end

end