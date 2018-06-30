function solMap = Map(solFrom, limits, orderFrom, orderTo, useFixUp, simplexType)

% Fine-to-coarse mapping
if orderFrom > orderTo    
    
    solMap = fineToCoarse(solFrom, limits, orderFrom, orderTo, useFixUp, simplexType);
    
    % Coarse-to-fine mapping
elseif orderFrom < orderTo
    
    solMap = coarseToFine(solFrom, limits, orderFrom, orderTo, useFixUp, simplexType);
    
    % One-to-one mapping
else
   
    error('Error: Attempting to map across the same orders.');
    
end

end