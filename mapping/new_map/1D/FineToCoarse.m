% Fine-to-coarse mapping
function solCoarse = FineToCoarse(solFrom, orderFrom, orderTo, limits)

useFixUp = 1;

%% Retrieve limits
psiMin = limits.psiMin;
psiMax = limits.psiMax;

%% Initial quadrature information

% Fine quadrature information
orderFine = orderFrom;
quadFine  = GenQuad(orderFine);
solFine   = solFrom;

% Coarse quadrature information
orderCoarse = orderFine - 1;
quadCoarse  = GenQuad(orderCoarse);

%% Go through each refinement level
while orderCoarse >= orderTo
    
    % Gather coarse quadrature information
    weightCoarse = quadCoarse.weight;
    muCoarse     = quadCoarse.mu;
    
    % Gather fine quadrature information
    weightFine = quadFine.weight;
    muFine     = quadFine.mu;
    
    % Initialize storage
    solCoarse = cell(size(muCoarse, 2), 1);
    
    % Go over each coarse LDFE region
    for i = 1 : size(muCoarse, 2)
        
        %% Apply mapping algorithm
        
        % Initialize A-matrix
        aMat = zeros(2);
        
        % Form A-matrix
        
        % Go through each coarse direction
        for j = 1 : 2
            aMat(1, j) = weightCoarse{i}(j);
            aMat(2, j) = weightCoarse{i}(j) * muCoarse{i}(j);
        end
        
        % Initialize b-vector
        bVec = zeros(2, 1);
        
        % From b-vector
        
        % Go through each fine LDFE region
        for k = 2 * i - 1 : 2 * i
            
            % Go through each fine direction
            for l = 1 : 2
                bVec(1) = bVec(1) + weightFine{k}(l) * solFine{k}(l);
                bVec(2) = bVec(2) + weightFine{k}(l) * solFine{k}(l) * muFine{k}(l);
            end
        end
        
        % Coarse solution
        solCoarse{i} = aMat \ bVec;
        
        %% Apply fix-up algorithm
        
        if useFixUp == 1
            
            % Initialize values
            deltaPhiR = 0; % Amount of 0th moment to reallocate
            counter   = 1; % Counter for good ordinates
            
            % Go through each coarse direction
            for j = 1 : 2
                
                % Check if mapped solution is greater than max limit
                if solCoarse{i}(j) > psiMax
                    
                    % Amount of angular flux to reallocate
                    deltaPsi = solCoarse{i}(j) - psiMax;
                    
                    % Amount of zeroth moment to reallocate
                    deltaPhiR = deltaPhiR + weightCoarse{i}(j) * deltaPsi;
                    
                    % Push bad ordinate to maximum limit
                    solCoarse{i}(j) = psiMax;
                    
                elseif solCoarse{i}(j) < psiMin
                    
                    % Amount of angular flux to reallocate
                    deltaPsi = solCoarse{i}(j) - psiMin;
                    
                    % Amount of zeroth moment to reallocate
                    deltaPhiR = deltaPhiR + weightCoarse{i}(j) * deltaPsi;
                    
                    % Push bad ordinate to minimum limit
                    solCoarse{i}(j) = psiMin;
                    
                else
                    
                    % Calculate 0th moment allowed from min limit
                    goodOrd(counter).deltaPhiMinAllowed = ...
                        weightCoarse{i}(j) * (solCoarse{i}(j) - psiMin);
                    
                    % Calculate 0th moment allowed from max limit
                    goodOrd(counter).deltaPhiMaxAllowed = ...
                        weightCoarse{i}(j) * (psiMax - solCoarse{i}(j));
                    
                    % Calculate 0th moment
                    goodOrd(counter).zerothMoment = ...
                        weightCoarse{i}(j) * solCoarse{i}(j);
                    
                    % Store ordinate direction and weight
                    goodOrd(counter).mu = muCoarse{i}(j);
                    goodOrd(counter).weight = weightCoarse{i}(j);
                    
                    % Store ordinate location
                    goodOrd(counter).i = i;
                    goodOrd(counter).j = j;
                    
                    % Update counter
                    counter = counter + 1;
                    
                end
                
            end
            
            % Number of good ordinates
            numGoodOrd = counter - 1;
            
            % If there is one good ordinate
            if numGoodOrd == 1
                
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
                solCoarse{goodOrd(1).i}(goodOrd(1).j) = ...
                    (goodOrd(1).zerothMoment + deltaPhiR) / goodOrd(1).weight;
                
                % We assume there are no mapped regions with all bad ordinates
            elseif numGoodOrd == 0
                
                error('Fixup algorithm failed.');
                
            end
            
            clearvars goodOrd;            
            
        end
        
    end
    
    % Update coarse and fine information
    if orderCoarse > orderTo
        
        % Fine quadrature order
        orderFine = orderFine - 1;
        
        % Fine quadrature
        quadFine = GenQuad(orderFine);
        
        % Fine quadrature solution
        solFine = solCoarse;
        
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