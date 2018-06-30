function solMap = Map(solFrom, limits, orderFrom, orderTo, useFixUp, LPtype)

%% Fine-to-coarse mapping
if orderFrom > orderTo    
    
    solMap = fineToCoarse(solFrom, limits, orderFrom, orderTo, useFixUp, LPtype);
    
    %% Coarse-to-fine mapping
elseif orderFrom < orderTo
    
    solMap = coarseToFine(solFrom, limits, orderFrom, orderTo, useFixUp, LPtype);
    
    %% One-to-one mapping
else
   
    error('Error: Attempting to map across the same order.');
    
end

end