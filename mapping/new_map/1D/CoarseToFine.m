% Coarse-to-fine mapping
function solFine = CoarseToFine(solFrom, orderFrom, orderTo, limits)

useFixUp = 0;

% User-specified limits
psiMin = limits.psiMin;
psiMax = limits.psiMax;

% Coarse quadrature information
orderCoarse = orderFrom;
quadCoarse  = GenQuad(orderCoarse);
solCoarse   = solFrom;

% Fine quadrature information
orderFine = orderFrom + 1;
quadFine  = GenQuad(orderFine);

% While fine order less than quadrature order you are mapping to
while orderFine <= orderTo
    
    % Gather coarse quadrature information
    weightCoarse = quadCoarse.weight;
    muCoarse     = quadCoarse.mu;
    basisCoarse  = quadCoarse.basis;
    
    % Gather fine quadrature information
    weightFine = quadFine.weight;
    muFine     = quadFine.mu;
    
    % Initialize storage
    solFine = cell(size(muFine, 2), 1);
    
    % Go through each coarse LDFE region
    for i = 1 : size(muCoarse, 2)
        
        %% Apply mapping algorithm
        
        % Initialize A-matrix
        aMat = zeros(2);
        
        % Form A-matrix
        
        % Go through each coarse direction
        for j = 1 : 2
            
            % Go through each fine LDFE region
            for k = (i - 1) * 2 + 1 : (i - 1) * 2 + 2
                
                % Go through each fine direction
                for m = 1 : 2
                    
                    aMat(1, j) = aMat(1, j) + weightFine{k}(m) * ...
                        (basisCoarse{i}(1, j) + ...
                        basisCoarse{i}(2, j) * ...
                        muFine{k}(m));
                    
                    aMat(2, j) = aMat(2, j) + weightFine{k}(m) * ...
                        muFine{k}(m) * ...
                        (basisCoarse{i}(1, j) + ...
                        basisCoarse{i}(2, j) * ...
                        muFine{k}(m));
                    
                end
                
            end
            
        end
        
        % Initialize b-vector
        bVec = zeros(2, 1);
        
        % Form b-vector
        
        % Go over each coarse direction
        for j = 1 : 2
            bVec(1) = bVec(1) + weightCoarse{i}(j) * ...
                solCoarse{i}(j);
            bVec(2) = bVec(2) + weightCoarse{i}(j) * ...
                muCoarse{i}(j) * solCoarse{i}(j);
        end
        
        % Calculate psi-tilde values
        coarseSolTilde = aMat \ bVec;
        
        % Calculate mapped fine solution
        
        % Go over each fine LDFE region
        for j = (i - 1) * 2 + 1 : (i - 1) * 2 + 2
            
            % Go over each fine direction
            for k = 1 : 2
                
                % Sum over each coarse basis function
                solFine{j}(k) = 0;
                for m = 1 : 2
                    solFine{j}(k) = solFine{j}(k) + ...
                        coarseSolTilde(m) * ...
                        (basisCoarse{i}(1, m) + ...
                        basisCoarse{i}(2, m) * ...
                        muFine{j}(k));
                end
                
            end
            
        end
        
        %% Apply mapping fix-up algorithm
        if useFixUp == 1
            
            % Initialize values
            deltaPhiR = 0; % Amount of 0th moment to reallocate
            deltaCurR = 0; % Amount of 1st moment to reallocate
            counter   = 1; % Counter for good ordinates
            
            % Go over each fine LDFE region
            for j = (i - 1) * 2 + 1 : (i - 1) * 2 + 2
                
                % Go over each fine direction
                for k = 1 : 2
                    
                    % Check if mapped solution is greater than max limit
                    if solFine{j}(k) > psiMax
                        
                        % Amount of angular flux to reallocate
                        deltaPsi = solFine{j}(k) - psiMax;
                        
                        % Amount of 0th and 1st moment to reallocate
                        deltaPhiR = deltaPhiR + weightFine{j}(k) * deltaPsi;
                        
                        deltaCurR = ...
                            deltaCurR + weightFine{j}(k) * muFine{j}(k) * deltaPsi;
                        
                        % Set bad ordinate to maximum limit
                        solFine{j}(k) = psiMax;
                        
                        % Check if mapped solution is less than min limit
                    elseif solFine{j}(k) < psiMin
                        
                        % Amount of angular flux to reallocate
                        deltaPsi = solFine{j}(k) - psiMin;
                        
                        % Amount of 0th and 1st moment to reallocate
                        deltaPhiR = deltaPhiR + weightFine{j}(k) * deltaPsi;
                        
                        deltaCurR = ...
                            deltaCurR + weightFine{j}(k) * muFine{j}(k) * deltaPsi;
                        
                        % Set bad ordinate to minimum limit
                        solFine{j}(k) = psiMin;
                        
                        % Check if mapped solution is within limits
                    else
                        
                        % Calculate 0th moment allowed from min limit
                        goodOrd(counter).deltaPhiMinAllowed = ...
                            weightFine{j}(k) * (solFine{j}(k) - psiMin);
                        
                        % Calculate 0th moment allowed from max limit
                        goodOrd(counter).deltaPhiMaxAllowed = ...
                            weightFine{j}(k) * (psiMax - solFine{j}(k));
                        
                        % Calculate 0th moment
                        goodOrd(counter).zerothMoment = ...
                            weightFine{j}(k) * solFine{j}(k);
                        
                        % Store ordinate direction and weight
                        goodOrd(counter).mu = muFine{j}(k);
                        goodOrd(counter).weight = weightFine{j}(k);
                        
                        % Store ordinate location
                        goodOrd(counter).j = j;
                        goodOrd(counter).k = k;
                        
                        % Update counter
                        counter = counter + 1;
                        
                    end
                    
                end
                
            end
            
            % Number of good ordinates
            numGoodOrd = counter - 1;
            
            % If there are three in-range angular flux solutions
            if numGoodOrd == 3
                
                % If we need to add 0th moment to the good ordinates
                if deltaPhiR > 0
                    
                    % Sort good ordinates based on 0th moment allowed from max
                    % limit
                    [~, sortedInd] = sort([goodOrd.deltaPhiMaxAllowed]);
                    sortedInd = fliplr(sortedInd);
                    goodOrd = goodOrd(sortedInd);
                    
                    % Make sure total 0th moment allowed is suitable
                    totalPhiAllowed = ...
                        goodOrd(1).deltaPhiMaxAllowed + ...
                        goodOrd(2).deltaPhiMaxAllowed + ...
                        goodOrd(3).deltaPhiMaxAllowed;
                    
                    if totalPhiAllowed < deltaPhiR
                        error('Error 1: Fixup algorithm failed');
                    end
                    
                    % If we need to remove 0th moment from the good ordinates
                else
                    
                    % Sort good ordinates based on 0th moment allowed from min
                    % limit
                    [~, sortedInd] = sort([goodOrd.deltaPhiMinAllowed]);
                    sortedInd = fliplr(sortedInd);
                    goodOrd = goodOrd(sortedInd);
                    
                    % Make sure total 0th moment allowed is suitable
                    totalPhiAllowed = ...
                        goodOrd(1).deltaPhiMinAllowed + ...
                        goodOrd(2).deltaPhiMinAllowed + ...
                        goodOrd(3).deltaPhiMinAllowed;
                    
                    if totalPhiAllowed < abs(deltaPhiR)
                        error('Error 2: Fixup algorithm failed');
                    end
                    
                end
                
                % Calculate alpha
                alpha = ...
                    (goodOrd(1).zerothMoment +  ...
                    goodOrd(2).zerothMoment) / ...
                    (goodOrd(1).zerothMoment +  ...
                    goodOrd(2).zerothMoment +  ...
                    goodOrd(3).zerothMoment);
                
                % Calculate beta
                beta = ...
                    (deltaCurR - alpha * deltaPhiR * goodOrd(2).mu - ...
                    (1 - alpha) * deltaPhiR * goodOrd(3).mu) / ...
                    (alpha * deltaPhiR * (goodOrd(1).mu - goodOrd(2).mu));
                
                % Adjust psi in top group
                solFine{goodOrd(1).j}(goodOrd(1).k) = ...
                    (goodOrd(1).zerothMoment + beta * alpha * deltaPhiR) / ...
                    goodOrd(1).weight;
                
                solFine{goodOrd(2).j}(goodOrd(2).k) = ...
                    (goodOrd(2).zerothMoment + (1 - beta) * alpha * deltaPhiR) / ...
                    goodOrd(2).weight;
                
                % Adjust psi in bottom group
                solFine{goodOrd(3).j}(goodOrd(3).k) = ...
                    (goodOrd(3).zerothMoment + (1 - alpha) * deltaPhiR) / ...
                    goodOrd(3).weight;
                
                % Group both good ordinates into the top set (no alpha value)
            elseif numGoodOrd == 2
                
                % If we need to add 0th moment to the good ordinates
                if deltaPhiR > 0
                    
                    % Sort good ordinates based on 0th moment allowed from max
                    % limit
                    [~, sortedInd] = sort([goodOrd.deltaPhiMaxAllowed]);
                    sortedInd = fliplr(sortedInd);
                    goodOrd = goodOrd(sortedInd);
                    
                    % Make sure total 0th moment allowed is suitable
                    totalPhiAllowed = ...
                        goodOrd(1).deltaPhiMaxAllowed + ...
                        goodOrd(2).deltaPhiMaxAllowed;
                    
                    if totalPhiAllowed < deltaPhiR
                        error('Error 1: Fixup algorithm failed');
                    end
                    
                    % If we need to remove 0th moment from the good ordinates
                else
                    
                    % Sort good ordinates based on 0th moment allowed from min
                    % limit
                    [~, sortedInd] = sort([goodOrd.deltaPhiMinAllowed]);
                    sortedInd = fliplr(sortedInd);
                    goodOrd = goodOrd(sortedInd);
                    
                    % Make sure total 0th moment allowed is suitable
                    totalPhiAllowed = ...
                        goodOrd(1).deltaPhiMinAllowed + ...
                        goodOrd(2).deltaPhiMinAllowed;
                    
                    if totalPhiAllowed < abs(deltaPhiR)
                        error('Error 2: Fixup algorithm failed');
                    end
                    
                end
                
                % Calculate beta
                beta = ...
                    (deltaCurR - deltaPhiR * goodOrd(2).mu) / ...
                    (deltaPhiR * (goodOrd(1).mu - goodOrd(2).mu));
                
                % Adjust psi in top group
                solFine{goodOrd(1).j}(goodOrd(1).k) = ...
                    (goodOrd(1).zerothMoment + beta * deltaPhiR) / goodOrd(1).weight;
                
                solFine{goodOrd(2).j}(goodOrd(2).k) = ...
                    (goodOrd(2).zerothMoment + (1 - beta) * deltaPhiR) / goodOrd(2).weight;
                
                % We can only preserve either the 0th or 1st moment
            elseif numGoodOrd == 1
                
                % If we need to add 0th moment to the good ordinates
                if deltaPhiR > 0
                    
                    % Make sure total 0th moment allowed is suitable
                    if (goodOrd(1).deltaPhiMaxAllowed < deltaPhiR && ...
                            abs(deltaPhiR) > eps)
                        
                        error('Error 1: Fixup algorithm failed');
                        
                    end
                    
                    % If we need to remove 0th moment from the good ordinates
                else
                    
                    if (goodOrd(1).deltaPhiMinAllowed < abs(deltaPhiR) && ...
                            abs(deltaPhiR) > eps)
                        
                        error('Error 2: Fixup algorithm failed');
                        
                    end
                    
                end
                
                % Adjust psi to perserve 0th moment
                solFine{goodOrd(1).j}(goodOrd(1).k) = ...
                    (goodOrd(1).zerothMoment + deltaPhiR) / goodOrd(1).weight;
                
                
                % We assume there are no mapped regions with all bad ordinates
            elseif numGoodOrd == 0
                
                error('Fixup algorithm failed.');
                
            end
            
            clearvars goodOrd;
            
        end
        
    end
    
    % Update coarse and fine information
    if orderFine < orderTo
        
        % Coarse quadrature order
        orderCoarse = orderCoarse + 1;
        
        % Coarse quadrature
        quadCoarse = GenQuad(orderCoarse);
        
        % Coarse quadrature solution
        solCoarse = solFine;
        
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