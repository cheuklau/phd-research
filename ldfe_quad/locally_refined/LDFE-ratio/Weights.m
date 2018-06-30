% Solve weights
function weights = Weights(integrations, constants)

% Initialize storage
weights = zeros(1, 4);

% Go through each basis function
for i = 1 : 4
    % Add contribution from each term
    for j = 1 : 4
        weights(i) = weights(i) + constants(j, i) * integrations(j);
    end
end

end