% Fine-to-coarse mapping
function solCoarse = fineToCoarse(solFrom, orderFrom, orderTo, mapping)

% Fine quadrature information
orderFine = orderFrom;
quadFine  = GenQuad(orderFine);
solFine   = solFrom;

% Coarse quadrature information
orderCoarse = orderFine - 1;
quadCoarse  = GenQuad(orderCoarse);

% Map coarser one level at a time
while orderCoarse >= orderTo
    
    % Coarse quadrature information
    weightCoarse = quadCoarse.weight;
    muCoarse     = quadCoarse.mu;
    
    % Fine quadrature information
    weightFine = quadFine.weight;
    muFine     = quadFine.mu;
    basisFine  = quadFine.basis;
    
    % Initialize storage
    numRegions = size(muCoarse, 2);
    solCoarse  = cell(numRegions, 1);
    
    % Preserve 0th and 1st moments
    if mapping == 1 
        
        % Go over each coarse LDFE region
        for i = 1 : numRegions
            
            % A-matrix
            aMat = zeros(2);
            for j = 1 : 2
                aMat(1, j) = weightCoarse{i}(j);
                aMat(2, j) = weightCoarse{i}(j) * muCoarse{i}(j);
            end
            
            % b-vector
            bVec = zeros(2, 1);
            for k = 2 * i - 1 : 2 * i
                for l = 1 : 2
                    bVec(1) = bVec(1) + weightFine{k}(l) * solFine{k}(l);
                    bVec(2) = bVec(2) + weightFine{k}(l) * solFine{k}(l) * muFine{k}(l);
                end
            end
            
            % Coarse solution
            solCoarse{i} = aMat \ bVec;
            
        end
                
        % Preserve 1st moment using LDFE functions
    elseif mapping == 2
        
        % Go over each coarse LDFE region
        for i = 1 : numRegions
            
            % Go over each coarse quadrature direction
            for j = 1 : 2
                
                % Fine LDFE region
                fineRange = 2 * i - 2 + j;
                
                % Mapped solution
                solCoarse{i}(j) = 0;
                for k = 1 : 2
                    solCoarse{i}(j) = ...
                        solCoarse{i}(j) + ...
                        solFine{fineRange}(k) * ...
                        (basisFine{fineRange}(1, k) +...
                         basisFine{fineRange}(2, k) * muCoarse{i}(j));
                end
                
            end
            
            % Conservation factor numerator
            numerator   = 0;            
            for j = 2 * i - 1 : 2 * i
                for k = 1 : 2
                    numerator = numerator + weightFine{j}(k) * ...
                        solFine{j}(k) * muFine{j}(k);
                end
            end
            
            % Conservation factor denominator
            denominator = 0;
            for j = 1 : 2
                denominator = denominator + weightCoarse{i}(j) * ...
                    solCoarse{i}(j) * muCoarse{i}(j);
            end
            
            % Apply conservation factor
            if abs(denominator) ~= 0
                factor = numerator / denominator;
                solCoarse{i} = solCoarse{i} * factor;
            end
            
        end
        
        % Preserve 1st moment using CDFE functions
    elseif mapping == 3
        
        % Go over each coarse LDFE region
        for i = 1 : numRegions
            
            % Go over each coarse quadrature direction
            for j = 1 : 2
                
                % Fine LDFE region
                fineRange = 2 * i - 2 + j;
                
                % Solution numerator
                numerator = 0;
                for k = 1 : 2
                    numerator = numerator + solFine{fineRange}(k) * ...
                        muFine{fineRange}(k) * weightFine{fineRange}(k);
                end
                
                % Solution denominator
                denominator = muCoarse{i}(j) * weightCoarse{i}(j);
                
                % Mapped solution
                solCoarse{i}(j) = numerator / denominator;
                
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