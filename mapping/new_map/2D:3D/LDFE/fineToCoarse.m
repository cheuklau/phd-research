function solCoarse = fineToCoarse(solFrom, limits, orderFrom, orderTo, useFixUp, LPtype)

%% Limits
psiMin = limits.psiMin;
psiMax = limits.psiMax;

%% Quadrature information

% Fine
orderFine = orderFrom;
quadFine  = GenQuad(orderFine);
solFine   = solFrom;

% Coarse
orderCoarse = orderFine - 1;
quadCoarse  = GenQuad(orderCoarse);

%% TEST --- Total number of iterations
tot_iters  = 0;
tot_fixups = 0;

%% Go through each refinement level
while orderCoarse >= orderTo
    
    %% Sub-square properties
    
    % Number of fine sub-squares per face
    numSubSqFine    = size(quadFine, 1) / (3 * 4);
    numSubSqFineRow = numSubSqFine ^ 0.5;
    
    % Number of coarse sub-squares per face
    numSubSqCoarse = size(quadCoarse, 1) / (3 * 4);
    
    %% Initialize fine solution storage
    solCoarse = cell(3, 1);
    
    %% Go through each face
    for iFace = 1 : 3
        
        %% Initialize fine sub-square positions
        deltaX  = 0;
        deltaY  = 0;
        
        %% Go through each coarse sub-square
        for i = 1 : numSubSqCoarse
            
            %% Apply mapping algorithm
            
            % Coarse sub-square directions and weights
            muCoarse  = zeros(4, 1);
            etaCoarse = zeros(4, 1);
            xiCoarse  = zeros(4, 1);
            wtCoarse  = zeros(4, 1);
            for j = 1 : 4
                iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                muCoarse(j)  = quadCoarse{iDir}(1);
                etaCoarse(j) = quadCoarse{iDir}(2);
                xiCoarse(j)  = quadCoarse{iDir}(3);
                wtCoarse(j)  = quadCoarse{iDir}(4);
            end
            
            % Fine direction indices
            lineFineTmp(1) = 4 * ((iFace - 1) * numSubSqFine + deltaX + deltaY * numSubSqFineRow) + 1;
            lineFineTmp(2) = lineFineTmp(1) + 4;
            lineFineTmp(3) = lineFineTmp(1) + numSubSqFineRow * 4;
            lineFineTmp(4) = lineFineTmp(3) + 4;
            lineFine(1) = lineFineTmp(4);
            lineFine(2) = lineFineTmp(2);
            lineFine(3) = lineFineTmp(1);
            lineFine(4) = lineFineTmp(3);
            
            % Initialize A-matrix
            aMat = zeros(4);
            
            % Go through each coarse direction
            for j = 1 : 4
                aMat(1, j) = wtCoarse(j);
                aMat(2, j) = wtCoarse(j) * muCoarse(j);
                aMat(3, j) = wtCoarse(j) * etaCoarse(j);
                aMat(4, j) = wtCoarse(j) * xiCoarse(j);
            end
            
            % Initialize b-vector
            bvec = zeros(4, 1);
            
            % Go though each fine sub-square
            for j = 1 : 4
                
                % Go through each fine direction
                for k = 0 : 3
                    
                    % Fine direction index
                    lineUse = lineFine(j) + k;
                    
                    % Add contribution to b-vector
                    bvec(1) = bvec(1) + ...
                        quadFine{lineUse}(4) * solFine(lineUse);
                    bvec(2) = bvec(2) + ...
                        quadFine{lineUse}(4) * quadFine{lineUse}(1) * solFine(lineUse);
                    bvec(3) = bvec(3) + ...
                        quadFine{lineUse}(4) * quadFine{lineUse}(2) * solFine(lineUse);
                    bvec(4) = bvec(4) + ...
                        quadFine{lineUse}(4) * quadFine{lineUse}(3) * solFine(lineUse);
                    
                end
                
            end
            
            % Solve for mapped coarse values
            solCoarse{iFace}{i} = aMat \ bvec;
            
            % Update fine sub-square position
            deltaX = deltaX + 2;
            if deltaX == numSubSqFineRow
                deltaX = 0;
                deltaY = deltaY + 2;
            end
            
            %% Apply fix-up algorithm
            
            if useFixUp == 1
                
                %% Calculate 0th and 1st moments from fine sub-squares
                
                % Initialize storage
                phiInd = zeros(4, 1);
                muInd  = zeros(4, 1);
                etaInd = zeros(4, 1);
                xiInd  = zeros(4, 1);
                
                % Go through each fine sub-square
                for j = 1 : 4
                    
                    % Go through each fine direction
                    for k = 0 : 3
                        
                        % Fine quadrature index
                        lineUse = lineFine(j) + k;
                        
                        % 0th moment for each fine direction
                        phiInd(j) = phiInd(j) + quadFine{lineUse}(4) * solFine(lineUse);
                        
                        % 1st moment for each fine direction
                        muInd(j)  = muInd(j)  + quadFine{lineUse}(4) * quadFine{lineUse}(1) * solFine(lineUse);
                        etaInd(j) = etaInd(j) + quadFine{lineUse}(4) * quadFine{lineUse}(2) * solFine(lineUse);
                        xiInd(j)  = xiInd(j)  + quadFine{lineUse}(4) * quadFine{lineUse}(3) * solFine(lineUse);
                        
                    end
                    
                end
                
                % Total 0th moment
                totPhi = sum(phiInd);
                
                % Total 1st moment
                totMu  = sum(muInd);
                totEta = sum(etaInd);
                totXi  = sum(xiInd);
                
                %% Check mapped solutions for extrema
                
                % Initialize out-of-range ordinate counter
                numBadOrd  = 0;
                
                % Go over each coarse direction
                for j = 1 : 4
                    
                    % If coarse solution is outside limits
                    if solCoarse{iFace}{i}(j) > psiMax || ...
                            solCoarse{iFace}{i}(j) < psiMin
                        
                        % Increment number of bad directions
                        numBadOrd = numBadOrd + 1;
                        
                    end
                    
                end
                
                %% Solve multi-objective problem
                if numBadOrd > 0
                    
                    %% A-matrix of constraint coefficients
                    
                    % Initialize A-matrix
                    Amat = zeros(12, 26);
                    
                    % 1st moment objectives and 0th moment constraint
                    for j = 1 : 4
                        
                        Amat(1, j) = wtCoarse(j) * muCoarse(j);
                        Amat(2, j) = wtCoarse(j) * etaCoarse(j);
                        Amat(3, j) = wtCoarse(j) * xiCoarse(j);
                        Amat(4, j) = wtCoarse(j);
                        
                    end
                    
                    % Main diagonal for slack and artifical variables
                    for j = 1 : 12
                        
                        Amat(j, 14 + j) = 1;
                        
                    end
                    
                    % Off diagonals for psi limits
                    for j = 1 : 4
                        
                        Amat(4 + j, j) = 1;
                        Amat(8 + j, j) = 1;
                        
                    end
                    
                    % Off diagonals for 1st moment deltas
                    for j = 1 : 3
                        
                        Amat(j, 4 + j) = -1;
                        Amat(j, 7 + j) = 1;
                        
                    end
                    
                    % Surplus variables for psi min limit
                    for j = 1 : 4
                        
                        Amat(4 + j, 10 + j) = -1;
                        
                    end
                    
                    %% b-vector of constraint right-hand sides
                    bVec = zeros(12, 1);
                    bVec(1) = totMu;
                    bVec(2) = totEta;
                    bVec(3) = totXi;
                    bVec(4) = totPhi;
                    for j = 5 : 8
                        bVec(j) = psiMin;
                    end
                    for j = 9 : 12
                        bVec(j) = psiMax;
                    end
                    
                    %% c-vector of cost function coefficients
                    cVec = zeros(1, 26);
                    for j = 5 : 10
                        cVec(j) = 1;
                    end
                    for j = 15 : 22
                        cVec(j) = 10;
                    end
                    
                    %% Initial feasible solution for Simplex methods
                    if LPtype ~= 4                        
                        bfs = zeros(1, 26);
                        bfs(15) = totMu;
                        bfs(16) = totEta;
                        bfs(17) = totXi;
                        bfs(18) = totPhi;
                        for j = 19 : 22
                            bfs(j) = psiMin;
                        end
                        for j = 23 : 26
                            bfs(j) = psiMax;
                        end                        
                    end
                    
                    %% Remove entries corresponding to zero right-hand side constraints
                    
                    row_to_delete = [];
                    col_to_delete = [];
                                       
                    for j = 1 : 12
                        
                        if bVec(j) == 0
                            
                            row_to_delete = [row_to_delete j];
                                                        
                            % 1st moment deltas and artificial
                            if j < 4
                                
                                col_to_delete = [col_to_delete 4+j 7+j];
                                
                                % 0th moment artificial
                            elseif j == 4
                                
                                col_to_delete = [col_to_delete 18];
                                
                                % Min limit surplus and artificial
                            elseif j < 9
                                
                                tmp = j - 4;
                                col_to_delete = [col_to_delete 10+tmp 18+tmp];
                                
                                % Max limit slack
                            else
                                
                                tmp = j - 8;
                                col_to_delete = [col_to_delete 22+tmp];
                                
                            end
                            
                        end
                        
                    end
                    Amat(row_to_delete,:)  = [];
                    Amat(:, col_to_delete) = [];
                    bVec(row_to_delete)    = [];
                    cVec(col_to_delete)    = [];
                    if LPtype ~= 4
                        bfs(col_to_delete) = [];
                    end
                                        
                    %% Call LP solver
                    
                    % Revised Simplex
                    if LPtype == 1
                        
                        [num_iters, x] = rsm(cVec, Amat, bVec, 1e-6, bfs);
                        
                        % Bartels-Golub
                    elseif LPtype == 2
                        
                        [num_iters, x] = bg(cVec, Amat, bVec, 1e-8, bfs);
                        
                        % Forrest-Tomlin
                    elseif LPtype == 3
                        
                        [num_iters, x] = ft(cVec, Amat, bVec, 1e-8, 1e-8, 1e-8, bfs);
                        
                        % Interior Point
                    elseif LPtype == 4
                        
                        [~, num_iters, x] = intpt(Amat, bVec, cVec');
                        
                    end
                    
                    %% Store results into mapped solution
                    
                    % Go through each coarse direction
                    for j = 1 : 4
                        
                        % Set near zero solutions to zero
                        if abs(x(j)) < 1e-10
                            
                            solCoarse{iFace}{i}(j) = 0;
                            
                            % Store remaining solutions
                        else
                            
                            solCoarse{iFace}{i}(j) = x(j);
                            
                        end
                        
                    end
                    
                    %% Test -- Update total fixups and iterations
                    tot_fixups = tot_fixups + 1;
                    tot_iters = tot_iters + num_iters;
                end
                
            end
            
        end
        
    end
    
    %% Update coarse and fine information
    if orderCoarse > orderTo
        
        % Next fine quadrature order
        orderFine = orderFine - 1;
        
        % Generate next fine quadrature
        quadFine = GenQuad(orderFine);
        
        % Next fine quadrature solution is current coarse solution
        numSubSqFine = size(quadFine, 1) / (3 * 4);
        counter = 1;
        for i = 1 : 3
            for j = 1 : numSubSqFine
                for k = 1 : 4
                    solFine(counter) = solCoarse{i}{j}(k);
                    counter = counter + 1;
                end
            end
        end
        
        % Next coarse quadrature order
        orderCoarse = orderCoarse - 1;
        
        % Generate next coarse quadrature
        quadCoarse = GenQuad(orderCoarse);
        
    else
        
        % Next coarse quadrature order
        orderCoarse = orderCoarse - 1;
        
    end
    
end

fprintf('Total number of fix-ups: %i \n', tot_fixups);
fprintf('Average number of iterations per fix-up: %e \n', tot_iters / tot_fixups);

end