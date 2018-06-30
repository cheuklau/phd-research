% Coarse-to-fine mapping
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
    
    % Gather coarse quadrature information
    weightCoarse = quadCoarse.weight;
    muCoarse     = quadCoarse.mu;
    basisCoarse  = quadCoarse.basis;
    
    % Gather fine quadrature information
    weightFine = quadFine.weight;
    muFine     = quadFine.mu;
    
    % Initialize storage
    numRegions = size(muFine, 2);
    solFine    = cell(numRegions, 1);
    
    % Preserve 0th and 1st moments
    if mapping == 1
        
        % Coarse LDFE region and direction trackers
        coarseRegion = 1;
        coarseDir   = 1;
        
        % Go over each fine LDFE region
        for i = 1 : numRegions
            
            % A-matrix
            aMat = zeros(2);
            for j = 1 : 2
                aMat(1, j) = weightFine{i}(j);
                aMat(2, j) = weightFine{i}(j) * muFine{i}(j);
            end
            
            % b-vector
            bVec = zeros(2, 1);
            bVec(1) = ...
                weightCoarse{coarseRegion}(coarseDir) * ...
                solCoarse{coarseRegion}(coarseDir);
            bVec(2) = ...
                weightCoarse{coarseRegion}(coarseDir) * ...
                solCoarse{coarseRegion}(coarseDir) * ...
                muCoarse{coarseRegion}(coarseDir);
            
            % Fine solution
            solFine{i} = aMat \ bVec;
            
            % Update coarse LDFE region and direction
            if mod(i, 2) == 0
                coarseRegion = coarseRegion + 1;
                coarseDir = 1;
            else
                coarseDir = coarseDir + 1;
            end
            
        end
        
        % Preserve 1st moment using LDFE functions
    elseif mapping == 2                
        
        % Coarse LDFE region trackers
        coarseRegion = 1;
        
        % Go over each fine LDFE region
        for i = 1 : numRegions
            
            % Go over each fine quadrature direction
            for j = 1 : 2
                
                % Mapped solution
                solFine{i}(j) = 0;
                for k = 1 : 2
                    solFine{i}(j) = ...
                        solFine{i}(j) + ...
                        solCoarse{coarseRegion}(k) * ...
                        (basisCoarse{coarseRegion}(1, k) + ...
                        basisCoarse{coarseRegion}(2, k) * muFine{i}(j));
                end
                
            end
            
            % Conservation factor
            if mod(i, 2) == 0
                
                % Conservation factor numerator
                numerator = 0;
                for j = 1 : 2
                    numerator = ...
                        numerator + ...
                        weightCoarse{coarseRegion}(j) * ...
                        solCoarse{coarseRegion}(j) * ...
                        muCoarse{coarseRegion}(j);
                end
                
                % Conservation factor denominator
                denominator = 0;
                for j = i - 1 : i
                    for k = 1 : 2
                        denominator = ...
                            denominator + ...
                            weightFine{j}(k) * ...
                            solFine{j}(k) * ...
                            muFine{j}(k);
                    end
                end
                
                % Apply conservation factor
                if abs(denominator) ~= 0
                    factor = numerator / denominator;
                    for j = i - 1 : i
                        solFine{j} = solFine{j} * factor;
                    end
                end
                
                % Update coarse LDFE region tracker
                coarseRegion = coarseRegion + 1;
                
            end
            
        end
        
        % Preserve 1st moment using CDFE functions
    elseif mapping == 3
        
        % Coarse LDFE region and direction trackers
        coarseRegion = 1;
        coarseDir    = 1;
        
        % Go over each fine LDFE region
        for i = 1 : numRegions
            
            % Solution numerator
            numerator = ...
                solCoarse{coarseRegion}(coarseDir) * ...
                muCoarse{coarseRegion}(coarseDir) * ...
                weightCoarse{coarseRegion}(coarseDir);
            
            % Solution denominator
            denominator = 0;
            for j = 1 : 2
                denominator = ...
                    denominator + ...
                    muFine{i}(j) * weightFine{i}(j);
            end
            
            % Mapped solution
            for j = 1 : 2
                solFine{i}(j) = numerator / denominator;
            end
            
            % Update coarse LDFE region and direction
            if mod(i, 2) == 0
                coarseRegion = coarseRegion + 1;
                coarseDir = 1;
            else
                coarseDir = coarseDir + 1;
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