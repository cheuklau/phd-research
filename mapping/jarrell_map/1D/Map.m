function solMap = Map(solFrom, orderFrom, orderTo, mapping)

% Fine-to-coarse mapping
if orderFrom > orderTo    
    
    solMap = fineToCoarse(solFrom, orderFrom, orderTo, mapping);
    
    % Coarse-to-fine mapping
elseif orderFrom < orderTo
    
    solMap = coarseToFine(solFrom, orderFrom, orderTo, mapping);
    
    % One-to-one mapping
else
   
    error('Error: Attempting to map across the same orders.');
    
end

end