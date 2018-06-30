function solMap = Map(solFrom, orderFrom, orderTo, limits)

% Fine-to-coarse mapping
if orderFrom > orderTo    
    
    solMap = FineToCoarse(solFrom, orderFrom, orderTo, limits);
    
    % Coarse-to-fine mapping
elseif orderFrom < orderTo
    
    solMap = CoarseToFine(solFrom, orderFrom, orderTo, limits);
    
    % One-to-one mapping
else
   
    error('Error: Attempting to map across the same orders.');
    
end

end