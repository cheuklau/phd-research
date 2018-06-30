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

%% TEST --- Total number of fix-ups and iterations
tot_iters  = 0;
tot_fixups = 0;

%% Go through each refinement level
while orderCoarse >= orderTo
    
    %% Sub-square properties
    
    % Number of fine sub-squares per face
    numSubSqFine    = size(quadFine, 1) / (3 * 9);
    numSubSqFineRow = numSubSqFine ^ 0.5;
    
    % Number of coarse sub-squares per face
    numSubSqCoarse = size(quadCoarse, 1) / (3 * 9);
    
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
            muCoarse  = zeros(9, 1);
            etaCoarse = zeros(9, 1);
            xiCoarse  = zeros(9, 1);
            wtCoarse  = zeros(9, 1);
            for j = 1 : 9
                iDir = (iFace - 1) * numSubSqCoarse * 9 + (i - 1) * 9 + j;
                muCoarse(j)  = quadCoarse{iDir}(1);
                etaCoarse(j) = quadCoarse{iDir}(2);
                xiCoarse(j)  = quadCoarse{iDir}(3);
                wtCoarse(j)  = quadCoarse{iDir}(4);
            end
            
            % Fine direction indices
            lineFine(1) = 9 * ((iFace - 1) * numSubSqFine + ...
                deltaX + deltaY * numSubSqFineRow) + 1;
            lineFine(2) = lineFine(1) + 9;
            lineFine(3) = lineFine(2) + 9;
            lineFine(4) = lineFine(1) + numSubSqFineRow * 9;
            lineFine(5) = lineFine(4) + 9;
            lineFine(6) = lineFine(5) + 9;
            lineFine(7) = lineFine(4) + numSubSqFineRow * 9;
            lineFine(8) = lineFine(7) + 9;
            lineFine(9) = lineFine(8) + 9;
            
            % Initialize A-matrix
            aMat = zeros(9);
            
            % Go through each coarse direction
            for j = 1 : 9
                mu = muCoarse(j);
                eta = etaCoarse(j);
                xi = xiCoarse(j);
                wt = wtCoarse(j);                
                aMat(1, j) = wt;
                aMat(2, j) = wt * mu;
                aMat(3, j) = wt * eta;
                aMat(4, j) = wt * xi;
                aMat(5, j) = wt * (mu ^ 2 - eta ^ 2);
                aMat(6, j) = wt * xi ^ 2;
                aMat(7, j) = wt * mu * eta;
                aMat(8, j) = wt * mu * xi;
                aMat(9, j) = wt * eta * xi;
            end
            
            % Initialize b-vector
            bvec = zeros(9, 1);
            
            % Go though each fine sub-square
            for j = 1 : 9
                
                % Go through each fine direction
                for k = 0 : 8
                    
                    % Fine direction index
                    lineUse = lineFine(j) + k;
                    
                    mu = quadFine{lineUse}(1);
                    eta = quadFine{lineUse}(2);
                    xi = quadFine{lineUse}(3);
                    wt = quadFine{lineUse}(4);
                    sol = solFine(lineUse);
                    
                    % Add contribution to b-vector
                    bvec(1) = bvec(1) + wt * sol;
                    bvec(2) = bvec(2) + wt * mu * sol;
                    bvec(3) = bvec(3) + wt * eta * sol;
                    bvec(4) = bvec(4) + wt * xi * sol;
                    bvec(5) = bvec(5) + wt * (mu ^ 2 - eta ^ 2) * sol;
                    bvec(6) = bvec(6) + wt * xi ^ 2 * sol;
                    bvec(7) = bvec(7) + wt * mu * eta * sol;
                    bvec(8) = bvec(8) + wt * mu * xi * sol;
                    bvec(9) = bvec(9) + wt * eta * xi * sol;
                    
                end
                
            end
            
            % Solve for mapped coarse values
            solCoarse{iFace}{i} = aMat \ bvec;
            
            % Update fine sub-square position
            deltaX = deltaX + 3;
            if deltaX == numSubSqFineRow
                deltaX = 0;
                deltaY = deltaY + 3;
            end
            
            %% Apply fix-up algorithm
            
            if useFixUp == 1
                
                %% Calculate 0th and 1st moments from fine sub-squares
                
                % Initialize storage
                phiInd     = zeros(9, 1);
                muInd      = zeros(9, 1);
                etaInd     = zeros(9, 1);
                xiInd      = zeros(9, 1);
                mu2eta2Ind = zeros(9, 1);
                xi2Ind     = zeros(9, 1);
                muetaInd   = zeros(9, 1);
                muxiInd    = zeros(9, 1);
                etaxiInd   = zeros(9, 1);
                
                % Go through each fine sub-square
                for j = 1 : 9
                    
                    % Go through each fine direction
                    for k = 0 : 8
                        
                        % Fine quadrature index
                        lineUse = lineFine(j) + k;
                        
                        mu  = quadFine{lineUse}(1);
                        eta = quadFine{lineUse}(2);
                        xi  = quadFine{lineUse}(3);
                        wt  = quadFine{lineUse}(4);
                        sol = solFine(lineUse);
                        
                        % 0th moment for each fine direction
                        phiInd(j) = phiInd(j) + wt * sol;
                        
                        % 1st moment for each fine direction
                        muInd(j)  = muInd(j)  + wt * mu * sol;
                        etaInd(j) = etaInd(j) + wt * eta * sol;
                        xiInd(j)  = xiInd(j)  + wt * xi * sol;                        
                        
                        % 2nd moment for each fine direction
                        mu2eta2Ind(j) = mu2eta2Ind(j) + wt * (mu ^ 2 - eta ^ 2) * sol;
                        xi2Ind(j)     = xi2Ind(j) + wt * xi ^ 2 * sol;
                        muetaInd(j)   = muetaInd(j) + wt * mu * eta * sol;
                        muxiInd(j)    = muxiInd(j) + wt * mu * xi * sol;
                        etaxiInd(j)   = etaxiInd(j) + wt * eta * xi * sol;
                        
                    end
                    
                end
                
                % Total 0th moment
                totPhi = sum(phiInd);
                
                % Total 1st moment
                totMu  = sum(muInd);
                totEta = sum(etaInd);
                totXi  = sum(xiInd);
                
                % Total 2nd moment
                totmu2eta2 = sum(mu2eta2Ind);
                totxi2 = sum(xi2Ind);
                totmueta = sum(muetaInd);
                totmuxi = sum(muxiInd);
                totetaxi = sum(etaxiInd);
                
                %% Check mapped solutions for extrema
                
                % Initialize out-of-range ordinate counter
                numBadOrd  = 0;
                
                % Go over each coarse direction
                for j = 1 : 9
                    
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
                    Amat = zeros(27, 61);
                    
                    % 0th moment constraint and 1st/2nd moment objectives
                    for j = 1 : 9
                    
                        mu  = muCoarse(j);
                        eta = etaCoarse(j);
                        xi  = xiCoarse(j);
                        wt  = wtCoarse(j);
                        
                        % 1st/2nd moment objectives
                        Amat(1, j) = wt * mu * eta;
                        Amat(2, j) = wt * mu * xi;
                        Amat(3, j) = wt * eta * xi;
                        Amat(4, j) = wt * xi ^ 2;
                        Amat(5, j) = wt * (mu ^ 2 - eta ^ 2);
                        Amat(6, j) = wt * mu;
                        Amat(7, j) = wt * eta;
                        Amat(8, j) = wt * xi;
                        
                        % 0th moment constraints                        
                        Amat(9, j) = wt;
                        
                    end
                    
                    % Main diagonal for artificial and slack variables
                    for j = 1 : 27
                       
                        Amat(j, 34 + j) = 1;
                        
                    end
                    
                    % Off diagonals for psi limits
                    for j = 1 : 9
                       
                        Amat(9 + j, j) = 1;
                        Amat(18 + j, j) = 1;
                        
                    end
                    
                    % Off diagonals for 1st/2nd moment deltas
                    for j = 1 : 8
                       
                        Amat(j, 9 + j) = -1;
                        Amat(j, 17 + j) = 1;
                        
                    end
                    
                    % Surplus variables for psi min limit
                    for j = 1 : 9
                       
                        Amat(9 + j, 25 + j) = -1;
                        
                    end
                    
                    %% b-vector of constraint right-hand sides
                    bVec = zeros(27, 1);
                    bVec(1) = totmueta;
                    bVec(2) = totmuxi;
                    bVec(3) = totetaxi;
                    bVec(4) = totxi2;
                    bVec(5) = totmu2eta2;
                    bVec(6) = totMu;
                    bVec(7) = totEta;
                    bVec(8) = totXi;
                    bVec(9) = totPhi;
                    for j = 10 : 18
                        bVec(j) = psiMin;
                    end
                    for j = 19 : 27
                        bVec(j) = psiMax;
                    end
                    
                    %% c-vector of cost function coefficients
                    cVec = zeros(1, 61);
                    for j = 10 : 25
                       cVec(j) = 1; 
                    end
                    for j = 35 : 52
                        cVec(j) = 10;
                    end
                    
                    %% Initial feasible solution for Simplex methods
                    if LPtype ~= 4
                       bfs = zeros(1, 61);
                       bfs(35) = totmueta;
                       bfs(36) = totmuxi;
                       bfs(37) = totetaxi;
                       bfs(38) = totxi2;
                       bfs(39) = totmu2eta2;
                       bfs(40) = totMu;
                       bfs(41) = totEta;
                       bfs(42) = totXi;
                       bfs(43) = totPhi;
                       for j = 44 : 52
                           bfs(j) = psiMin;
                       end
                       for j = 53 : 61
                           bfs(j) = psiMax;
                       end
                    end
                    
                    %% Remove entries corresponding to zero right-hand side constraints
                    row_to_delete = [];
                    col_to_delete = [];
                                       
                    for j = 1 : 27
                        
                        if bVec(j) == 0
                            
                            row_to_delete = [row_to_delete j];
                                                        
                            % 1st/2nd moment deltas and artificial
                            if j < 9
                                
                                col_to_delete = [col_to_delete 9+j 17+j 34+j];
                                
                                % 0th moment artificial
                            elseif j == 9
                                
                                col_to_delete = [col_to_delete 43];
                                
                                % Min limit surplus and artificial
                            elseif j < 19
                                
                                tmp = j - 9;
                                col_to_delete = [col_to_delete 25+tmp 43+tmp];
                                
                                % Max limit slack
                            else
                                
                                tmp = j - 18;
                                col_to_delete = [col_to_delete 52+tmp];
                                
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
                    for j = 1 : 9
                        
                        % Set near zero solutions to zero
                        if abs(x(j)) < 1e-10
                            
                            solCoarse{iFace}{i}(j) = 0;
                            
                            % Store remaining solutions
                        else
                            
                            solCoarse{iFace}{i}(j) = x(j);
                            
                        end
                        
                    end
                    
                    %% Test -- Update total fix-ups and iterations
                    tot_fixups = tot_fixups + 1;
                    tot_iters  = tot_iters + num_iters;
                    
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
        numSubSqFine = size(quadFine, 1) / (3 * 9);
        counter = 1;
        for i = 1 : 3
            for j = 1 : numSubSqFine
                for k = 1 : 9
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