function solFine = coarseToFine(solFrom, limits, orderFrom, orderTo, useFixUp, LPtype)

%% Limits
psiMin = limits.psiMin;
psiMax = limits.psiMax;

%% Quadrature information

% Coarse
orderCoarse = orderFrom;
quadCoarse  = GenQuad(orderCoarse);
solCoarse   = solFrom;

% Fine
orderFine = orderFrom + 1;
quadFine  = GenQuad(orderFine);

%% TEST --- Total number of iterations and fix-ups
tot_iters  = 0;
tot_fixups = 0;

%% Go through each refinement level
while orderFine <= orderTo
    
    %% Sub-square properties
    
    % Number of fine sub-squares per face
    numSubSqFine    = size(quadFine, 1) / (3 * 4);
    numSubSqFineRow = numSubSqFine ^ 0.5;
    
    % Number of coarse sub-squares per face
    numSubSqCoarse = size(quadCoarse, 1) / (3 * 4);
    
    %% Initialize fine solution storage
    solFine = cell(3, 1);
    
    %% Go through each face
    for iFace = 1 : 3
        
        %% Initialize fine sub-square positions
        deltaX = 0;
        deltaY = 0;
        
        %% Go through each coarse sub-square
        for i = 1 : numSubSqCoarse
            
            %% Apply mapping algorithm
            
            % Fine sub-square indices
            subSqFine(3) = deltaX + deltaY * numSubSqFineRow + 1;
            subSqFine(2) = subSqFine(3) + 1;
            subSqFine(4) = subSqFine(3) + numSubSqFineRow;
            subSqFine(1) = subSqFine(4) + 1;
            
            % Fine direction indices
            lineFine(3) = 4 * ((iFace - 1) * numSubSqFine + deltaX + deltaY * numSubSqFineRow) + 1;
            lineFine(2) = lineFine(3) + 4;
            lineFine(4) = lineFine(3) + numSubSqFineRow * 4;
            lineFine(1) = lineFine(4) + 4;
            
            % Initialize A-matrix
            aMat = zeros(4);
            
            % Go through each coarse direction
            for j = 1 : 4
                
                % Coarse direction index
                iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                
                % Go through each corresponding fine sub-square
                for k = 1 : 4
                    
                    % Go through each fine direction
                    for m = 1 : 4
                        
                        % Current fine sub-square data line
                        lineUse = lineFine(k) + m - 1;
                        
                        % Add contribution to A-matrix
                        aMat(1, j) = aMat(1, j) + ...
                            quadFine{lineUse}(4) * ...
                            (quadCoarse{iDir}(5) + ...
                            quadCoarse{iDir}(6) * quadFine{lineUse}(1) + ...
                            quadCoarse{iDir}(7) * quadFine{lineUse}(2) + ...
                            quadCoarse{iDir}(8) * quadFine{lineUse}(3));
                        aMat(2, j) = aMat(2, j) + ...
                            quadFine{lineUse}(4) * quadFine{lineUse}(1) * ...
                            (quadCoarse{iDir}(5) + ...
                            quadCoarse{iDir}(6) * quadFine{lineUse}(1) + ...
                            quadCoarse{iDir}(7) * quadFine{lineUse}(2) + ...
                            quadCoarse{iDir}(8) * quadFine{lineUse}(3));
                        aMat(3, j) = aMat(3, j) + ...
                            quadFine{lineUse}(4) * quadFine{lineUse}(2) * ...
                            (quadCoarse{iDir}(5) + ...
                            quadCoarse{iDir}(6) * quadFine{lineUse}(1) + ...
                            quadCoarse{iDir}(7) * quadFine{lineUse}(2) + ...
                            quadCoarse{iDir}(8) * quadFine{lineUse}(3));
                        aMat(4, j) = aMat(4, j) + ...
                            quadFine{lineUse}(4) * quadFine{lineUse}(3) * ...
                            (quadCoarse{iDir}(5) + ...
                            quadCoarse{iDir}(6) * quadFine{lineUse}(1) + ...
                            quadCoarse{iDir}(7) * quadFine{lineUse}(2) + ...
                            quadCoarse{iDir}(8) * quadFine{lineUse}(3));
                        
                    end
                    
                end
                
            end
            
            % Initialize b-vector
            bvec = zeros(4, 1);
            
            % Go through each coarse direction
            for j = 1 : 4
                
                % Coarse direction index
                iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                
                % Add contribution to b-vector
                bvec(1) = bvec(1) + ...
                    quadCoarse{iDir}(4) * solCoarse(iDir);
                bvec(2) = bvec(2) + ...
                    quadCoarse{iDir}(1) * quadCoarse{iDir}(4) * solCoarse(iDir);
                bvec(3) = bvec(3) + ...
                    quadCoarse{iDir}(2) * quadCoarse{iDir}(4) * solCoarse(iDir);
                bvec(4) = bvec(4) + ...
                    quadCoarse{iDir}(3) * quadCoarse{iDir}(4) * solCoarse(iDir);
                
            end
            
            % Solve for psi-tilde values
            coarseSolTilde = aMat \ bvec;
            
            % Go over each fine sub-square
            for j = 1 : 4
                
                % Go over each fine direction
                for k = 1 : 4
                    
                    % Initialize values
                    solFine{iFace}{subSqFine(j)}(k) = 0;
                    
                    % Fine quadrature index
                    lineUse = lineFine(j) + k - 1;
                    
                    % Go over each coarse direction
                    for m = 1 : 4
                        
                        % Coarse quadrature index
                        iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + m;
                        
                        % Mapped fine solution
                        solFine{iFace}{subSqFine(j)}(k) = ...
                            solFine{iFace}{subSqFine(j)}(k) + ...
                            coarseSolTilde(m) * ...
                            (quadCoarse{iDir}(5) + ...
                            quadCoarse{iDir}(6) * quadFine{lineUse}(1) + ...
                            quadCoarse{iDir}(7) * quadFine{lineUse}(2) + ...
                            quadCoarse{iDir}(8) * quadFine{lineUse}(3));
                        
                    end
                    
                end
                
            end
            
            % Update fine sub-square position
            deltaX = deltaX + 2;
            if deltaX == numSubSqFineRow
                deltaX = 0;
                deltaY = deltaY + 2;
            end
            
            %% Apply fix-up algorithm
            
            if useFixUp == 1
                
                %% Check mapped solution for extrema and calculate net change in scalar flux
                phi_net   = 0;
                phi_bad   = 0;
                mu_bad    = 0;
                eta_bad   = 0;
                xi_bad    = 0;
                numBadOrd = 0;
                good_ord = zeros(1, 16);
                
                for j = 1 : 4
                    
                    for k = 1 : 4
                        
                        lineUse = lineFine(j) + k - 1;
                        
                        if solFine{iFace}{subSqFine(j)}(k) > psiMax
                            
                            phi_net = phi_net + ...
                                quadFine{lineUse}(4) * (solFine{iFace}{subSqFine(j)}(k) - psiMax);
                            
                            solFine{iFace}{subSqFine(j)}(k) = psiMax;
                            
                            phi_bad = phi_bad + quadFine{lineUse}(4) * psiMax;
                            mu_bad  = mu_bad  + quadFine{lineUse}(4) * quadFine{lineUse}(1) * psiMax;
                            eta_bad = eta_bad + quadFine{lineUse}(4) * quadFine{lineUse}(2) * psiMax;
                            xi_bad  = xi_bad  + quadFine{lineUse}(4) * quadFine{lineUse}(3) * psiMax;
                            
                            numBadOrd = numBadOrd + 1;
                            
                            
                        elseif solFine{iFace}{subSqFine(j)}(k) < psiMin
                            
                            phi_net = phi_net + ...
                                quadFine{lineUse}(4) * (solFine{iFace}{subSqFine(j)}(k) - psiMin);
                            
                            solFine{iFace}{subSqFine(j)}(k) = psiMin;
                            
                            phi_bad = phi_bad + quadFine{lineUse}(4) * psiMin;
                            mu_bad  = mu_bad  + quadFine{lineUse}(4) * quadFine{lineUse}(1) * psiMin;
                            eta_bad = eta_bad + quadFine{lineUse}(4) * quadFine{lineUse}(2) * psiMin;
                            xi_bad  = xi_bad  + quadFine{lineUse}(4) * quadFine{lineUse}(3) * psiMin;
                            
                            numBadOrd = numBadOrd + 1;
                            
                        else
                            
                            good_ord((j - 1) * 4 + k) = 1;
                            
                        end
                        
                    end
                    
                end
                
                %% Solve multi-objective problem
                if numBadOrd > 0
                    
                    
                    %% Calculate 0th and 1st moments from coarse sub-square
                    
                    % Initialize storage
                    phi_total = 0;
                    mu_total  = 0;
                    eta_total = 0;
                    xi_total  = 0;
                    
                    % Go through each coarse direction
                    for j = 1 : 4
                        
                        % Coarse quadrature index
                        iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                        
                        % 0th moment for each coarse direction
                        phi_total = phi_total + quadCoarse{iDir}(4) * solCoarse(iDir);
                        
                        % 1st moment for each coarse direction
                        mu_total  = mu_total  + quadCoarse{iDir}(4) * quadCoarse{iDir}(1) * solCoarse(iDir);
                        eta_total = eta_total + quadCoarse{iDir}(4) * quadCoarse{iDir}(2) * solCoarse(iDir);
                        xi_total  = xi_total  + quadCoarse{iDir}(4) * quadCoarse{iDir}(3) * solCoarse(iDir);
                        
                    end
                    
                    %% Determine good ordinate limits for positive net scalar flux
                    if phi_net > 0
                        
                        % Initialize counter
                        counter = 1;
                        
                        % Reset temporary scalar flux
                        phi_temp = 0;
                        
                        psi_min_ind = zeros(1, 16);
                        
                        % Go through each fine sub-square
                        for j = 1 : 4
                            
                            % Go through each direction of current
                            % sub-square
                            for k = 1 : 4
                                
                                lineUse = lineFine(j) + k - 1;
                                
                                % If current direction is a good ordinate
                                if good_ord(counter) == 1
                                    
                                    % Set the lower limit as the original
                                    % mapped solution
                                    psi_min_ind(counter) = solFine{iFace}{subSqFine(j)}(k);
                                    
                                    
                                else
                                    
                                    % Dummy lower limit
                                    psi_min_ind(counter) = -1;
                                    
                                end
                                
                                % Increment scalar flux sum of current fine
                                % sub-square
                                phi_temp = phi_temp + ...
                                    quadFine{lineUse}(4) * solFine{iFace}{subSqFine(j)}(k);
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                        counter = 1;
                        
                        psi_max_ind = zeros(1, 16);
                        
                        % Go through each fine sub-square
                        for j = 1 : 4
                            
                            for k = 1 : 4
                                
                                if good_ord(counter) == 1
                                    
                                    psi_max_temp = (phi_total - (phi_temp - quadFine{lineUse}(4) * psi_min_ind(counter))) / ...
                                        quadFine{lineUse}(4);
                                    
                                    if psiMax < psi_max_temp
                                        
                                        psi_max_ind(counter) = psiMax;
                                        
                                    else
                                        
                                        psi_max_ind(counter) = psi_max_temp;
                                        
                                    end
                                    
                                else
                                    
                                    psi_max_ind(counter) = -1.0;
                                    
                                end
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                    elseif phi_net < 0
                        
                        % Initialize counter
                        counter = 1;
                        
                        % Reset temporary scalar flux
                        phi_temp = 0;
                        
                        psi_max_ind = zeros(1, 16);
                        
                        % Go through each fine sub-square
                        for j = 1 : 4
                            
                            % Go through each direction of current
                            % sub-square
                            for k = 1 : 4
                                
                                lineUse = lineFine(j) + k - 1;
                                
                                % If current direction is a good ordinate
                                if good_ord(counter) == 1
                                    
                                    % Set the upper limit as the original
                                    % mapped solution
                                    psi_max_ind(counter) = solFine{iFace}{subSqFine(j)}(k);
                                    
                                    
                                else
                                    
                                    % Dummy lower limit
                                    psi_max_ind(counter) = -1;
                                    
                                end
                                
                                % Increment scalar flux sum of current fine
                                % sub-square
                                phi_temp = phi_temp + ...
                                    quadFine{lineUse}(4) * solFine{iFace}{subSqFine(j)}(k);
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                        counter = 1;
                        
                        psi_min_ind = zeros(1, 16);
                        
                        % Go through each fine sub-square
                        for j = 1 : 4
                            
                            for k = 1 : 4
                                
                                if good_ord(counter) == 1
                                    
                                    psi_min_temp = (phi_total - (phi_temp - quadFine{lineUse}(4) * psi_max_ind(counter))) / ...
                                        quadFine{lineUse}(4);
                                    
                                    if psiMin > psi_min_temp
                                        
                                        psi_min_ind(counter) = psiMin;
                                        
                                    else
                                        
                                        psi_min_ind(counter) = psi_min_temp;
                                        
                                    end
                                    
                                else
                                    
                                    psi_min_ind(counter) = -1.0;
                                    
                                end
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                    end
                    
                    %% Initialize RSM variables
                    numGoodOrd = 16 - numBadOrd;
                    
                    %% A-matrix of constraint coefficients
                    
                    % Initialize A-matrix
                    Amat = zeros(4 + 2 * numGoodOrd, 10 + 4 * numGoodOrd);
                    
                    % 1st moment objectives
                    counter = 1;
                    
                    for j = 1 : 4
                        
                        for k = 1 : 4
                            
                            lineUse = lineFine(j) + k - 1;
                            
                            if good_ord((j - 1) * 4 + k) == 1
                                
                                Amat(1, counter) = quadFine{lineUse}(4) * quadFine{lineUse}(1);
                                Amat(2, counter) = quadFine{lineUse}(4) * quadFine{lineUse}(2);
                                Amat(3, counter) = quadFine{lineUse}(4) * quadFine{lineUse}(3);
                                Amat(4, counter) = quadFine{lineUse}(4);
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                    end
                    
                    % Main diagonal for slack and artificial variables
                    for j = 1 : 4 + 2 * numGoodOrd
                        
                        Amat(j, 2 * numGoodOrd + 6 + j) = 1;
                        
                    end
                    
                    % Off diagonals for psi limits
                    for j = 1 : numGoodOrd
                        
                        Amat(4 + j, j)  = 1;
                        Amat(4 + numGoodOrd + j, j) = 1;
                        
                    end
                    
                    % Off diagonals for 1st moment deltas
                    for j = 1 : 3
                        
                        Amat(j, numGoodOrd + j) = -1;
                        Amat(j, numGoodOrd + 3 + j) = 1;
                        
                    end
                    
                    % Surplus variables for psi min limit
                    for j = 1 : numGoodOrd
                        
                        Amat(4 + j, numGoodOrd + 6 + j) = -1;
                        
                    end
                    
                    %% b-vector of constraint right-hand sides
                    bVec = zeros(4 + 2 * numGoodOrd, 1);
                    bVec(1) = mu_total - mu_bad;
                    bVec(2) = eta_total - eta_bad;
                    bVec(3) = xi_total - xi_bad;
                    bVec(4) = phi_total - phi_bad;
                    
                    counter = 1;
                    
                    for j = 1 : 4
                        
                        for k = 1 : 4
                            
                            if good_ord((j - 1) * 4 + k) == 1
                                
                                bVec(4 + counter) = psi_min_ind((j - 1) * 4 + k);
                                
                                bVec(4 + numGoodOrd + counter) = psi_max_ind((j - 1) * 4 + k);
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                    end
                    
                    %% c-vector of cost function coefficients
                    cVec = zeros(1, 10 + 4 * numGoodOrd);
                    for j = numGoodOrd + 1 : numGoodOrd + 6
                        cVec(j) = 1;
                    end
                    for j = 7 + 2 * numGoodOrd : 10 + 3 * numGoodOrd
                        cVec(j) = 10;
                    end
                    
                    %% Initial feasible solution for Simplex method
                    if LPtype ~= 4
                        bfs = zeros(1, 10 + 4 * numGoodOrd);
                        bfs(2 * numGoodOrd + 7) = mu_total - mu_bad;
                        bfs(2 * numGoodOrd + 8) = eta_total - eta_bad;
                        bfs(2 * numGoodOrd + 9) = xi_total - xi_bad;
                        bfs(2 * numGoodOrd + 10) = phi_total - phi_bad;
                        
                        counter = 1;
                        
                        for j = 1 : 4
                            
                            for k = 1 : 4
                                
                                if(good_ord((j - 1) * 4 + k)) == 1
                                    
                                    bfs(2 * numGoodOrd + 10 + counter) = psi_min_ind((j - 1) * 4 + k);
                                    
                                    bfs(3 * numGoodOrd + 10 + counter) = psi_max_ind((j - 1) * 4 + k);
                                    
                                    counter = counter + 1;
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                    Amat_tmp = Amat;
                    cVec_tmp = cVec;
                    bVec_tmp = bVec;
                    bfs_tmp = bfs;
                    %% Remove entries corresponding to zero right-hand side constraints
                    
                    % Initialize rows and columns to delete
                    row_to_delete = [];
                    col_to_delete = [];
                    
                    % Go through each row
                    for j = 1 : 4 + 2 * numGoodOrd
                        
                        if bVec(j) == 0
                            
                            % Row to delete
                            row_to_delete = [row_to_delete j];
                            
                            % Column to delete
                            
                            % 1st moment delta and artificial variabls
                            if j <= 3
                                
                                col_to_delete = [col_to_delete j+numGoodOrd j+numGoodOrd+3 j+2*numGoodOrd+6];
                                
                                % 0th moment artificial variable
                            elseif j == 4
                                
                                col_to_delete = [col_to_delete 4+2*numGoodOrd+6];
                                
                                % Min limit surplus and artificial variables
                            elseif j <= numGoodOrd + 4
                                
                                temp = j - 4;
                                
                                col_to_delete = [col_to_delete temp+numGoodOrd+6 temp+2*numGoodOrd+6+4];
                                
                                % Max limit slack variable
                            else
                                
                                temp = j - numGoodOrd - 4;
                                
                                col_to_delete = [col_to_delete temp+2*numGoodOrd+6+4+numGoodOrd];
                                
                            end
                            
                        end
                        
                    end
                    
                    % Delete row and column of A-matrix
                    Amat(row_to_delete, :) = [];
                    Amat(:, col_to_delete) = [];
                    
                    % Delete b-vector elements
                    bVec(row_to_delete)    = [];
                    
                    % Delete c-vector elements
                    cVec(col_to_delete)    = [];
                    
                    % Delete initial feasible vector element for Simplex
                    if LPtype ~= 4
                        bfs(col_to_delete) = [];
                    end
                    
                    %% Call LP solver
                    
                    % Revised Simplex
                    if LPtype == 1
                        
                        % Test
                        B_indices = find(bfs);
                        if size(Amat(:,B_indices), 1) ~= size(Amat(:, B_indices), 2)
                            
                            disp('here');
                            
                        end
                        
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
                    counter = 1;
                    
                    for j = 1 : 4
                        
                        for k = 1 : 4
                            
                            if good_ord((j - 1) * 4 + k) == 1
                                
                                if abs(x(counter)) < 1e-10
                                    
                                    solFine{iFace}{subSqFine(j)}(k) = 0;
                                    
                                else
                                    
                                    solFine{iFace}{subSqFine(j)}(k) = x(counter);
                                    
                                end
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                    end
                    
                    %% Test -- Update total fix-ups and iterations
                    tot_fixups = tot_fixups + 1;
                    tot_iters = tot_iters + num_iters;
                    
                end
                
                
                
                %{
                %% Calculate 0th and 1st moments from coarse sub-square
                
                % Initialize storage
                phiInd = zeros(4, 1);
                muInd  = zeros(4, 1);
                etaInd = zeros(4, 1);
                xiInd  = zeros(4, 1);
                
                % Go through each coarse direction
                for j = 1 : 4
                    
                    % Coarse quadrature index
                    iDir = (iFace - 1) * numSubSqCoarse * 4 + (i - 1) * 4 + j;
                    
                    % 0th moment for each coarse direction
                    phiInd(j) = quadCoarse{iDir}(4) * solCoarse(iDir);
                    
                    % 1st moment for each coarse direction
                    muInd(j)  = quadCoarse{iDir}(4) * quadCoarse{iDir}(1) * solCoarse(iDir);
                    etaInd(j) = quadCoarse{iDir}(4) * quadCoarse{iDir}(2) * solCoarse(iDir);
                    xiInd(j)  = quadCoarse{iDir}(4) * quadCoarse{iDir}(3) * solCoarse(iDir);
                    
                end
                
                % Total 1st moment
                totMu  = sum(muInd);
                totEta = sum(etaInd);
                totXi  = sum(xiInd);
                
                %% Check mapped solutions for extrema
                
                % Initialize out-of-range ordinate counter
                numBadOrd  = 0;
                
                % Go over each fine sub-square
                for j = 1 : 4
                    
                    % Go over each fine direction
                    for k = 1 : 4
                        
                        % If mapped solution outside limits
                        if solFine{iFace}{subSqFine(j)}(k) > psiMax || ...
                                solFine{iFace}{subSqFine(j)}(k) < psiMin
                            
                            % Increment number of bad directions
                            numBadOrd = numBadOrd + 1;
                            
                        end
                        
                    end
                    
                end
                
                %% Solve multi-objective problem
                if numBadOrd > 0
                    
                    %% A-matrix of constraint coefficients
                    
                    % Initialize A-matrix
                    Amat = zeros(39, 77);
                    
                    % 1st moment objectives
                    counter = 1;
                    for j = 1 : 4
                        
                        for k = 1 : 4
                            
                            lineUse = lineFine(j) + k - 1;
                            
                            Amat(1, counter) = quadFine{lineUse}(4) * quadFine{lineUse}(1);
                            Amat(2, counter) = quadFine{lineUse}(4) * quadFine{lineUse}(2);
                            Amat(3, counter) = quadFine{lineUse}(4) * quadFine{lineUse}(3);
                            
                            counter = counter + 1;
                            
                        end
                        
                    end
                    
                    % 0th moment constraints
                    for j = 1 : 4
                        
                        for k = 1 : 4
                            
                            lineUse = lineFine(j) + k - 1;
                            
                            Amat(3 + j, (j - 1) * 4 + k) = quadFine{lineUse}(4);
                            
                        end
                        
                    end
                    
                    % Main diagonal for slack and artificial variables
                    for j = 1 : 39
                        
                        Amat(j, 38 + j) = 1;
                        
                    end
                    
                    % Off diagonals for psi limits
                    for j = 1 : 16
                        
                        Amat(7 + j, j)  = 1;
                        Amat(23 + j, j) = 1;
                        
                    end
                    
                    % Off diagonals for 1st moment deltas
                    for j = 1 : 3
                        
                        Amat(j, 16 + j) = -1;
                        Amat(j, 19 + j) = 1;
                        
                    end
                    
                    % Surplus variables for psi min limit
                    for j = 1 : 16
                        
                        Amat(7 + j, 22 + j) = -1;
                        
                    end
                    
                    %% b-vector of constraint right-hand sides
                    bVec = zeros(39, 1);
                    bVec(1) = totMu;
                    bVec(2) = totEta;
                    bVec(3) = totXi;
                    for j = 1 : 4
                        bVec(3 + j) = phiInd(j);
                    end
                    for j = 1 : 16
                        bVec(7 + j) = psiMin;
                        bVec(23 + j) = psiMax;
                    end
                    
                    %% c-vector of cost function coefficients
                    cVec = zeros(1, 77);
                    for j = 17 : 22
                        cVec(j) = 1;
                    end
                    for j = 39 : 61
                        cVec(j) = 10;
                    end
                    
                    %% Initial feasible solution for Simplex method
                    if LPtype ~= 4
                        bfs = zeros(1, 77);
                        bfs(39) = totMu;
                        bfs(40) = totEta;
                        bfs(41) = totXi;
                        for j = 1 : 4
                            bfs(41 + j) = phiInd(j);
                        end
                        for j = 46 : 61
                            bfs(j) = psiMin;
                        end
                        for j = 62 : 77
                            bfs(j) = psiMax;
                        end
                    end
                    
                    %% Remove entries corresponding to zero right-hand side constraints
                    
                    % Initialize rows and columns to delete
                    row_to_delete = [];
                    col_to_delete = [];
                    
                    % Go through each row
                    for j = 1 : 39
                        
                        if bVec(j) == 0
                            
                            % Row to delete
                            row_to_delete = [row_to_delete j];
                                           
                            % Column to delete
                            
                            % 1st moment delta and artificial variabls
                            if j < 4
                                
                                col_to_delete = [col_to_delete 16+j 19+j];
                                
                                % 0th moment artificial variable
                            elseif j < 8
                                
                                tmp = j - 3;
                                col_to_delete = [col_to_delete 41+tmp];
                                
                                % Min limit surplus and artificial variables
                            elseif j < 24
                                
                                tmp = j - 7;
                                col_to_delete = [col_to_delete 22+tmp 45+tmp];
                                
                                % Max limit slack variable
                            else
                                
                                tmp = j - 23;
                                col_to_delete = [col_to_delete 61+tmp];
                                
                            end
                            
                        end
                        
                    end
                    
                    % Delete row and column of A-matrix
                    Amat(row_to_delete, :) = [];
                    Amat(:, col_to_delete) = [];
                    
                    % Delete b-vector elements
                    bVec(row_to_delete)    = [];
                    
                    % Delete c-vector elements
                    cVec(col_to_delete)    = [];
                    
                    % Delete initial feasible vector element for Simplex
                    if LPtype ~= 4
                        bfs(col_to_delete) = [];
                    end
                    
                    %% Call LP solver
                    
                    % Revised Simplex
                    if LPtype == 1
                        
                        [num_iters, x] = rsm(cVec, Amat, bVec, 1e-10, bfs);
                        
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
                    counter = 1;
                    for j = 1 : 4
                        
                        for k = 1 : 4
                            
                            if abs(x(counter)) < 1e-10
                                
                                solFine{iFace}{subSqFine(j)}(k) = 0;
                                
                            else
                                
                                solFine{iFace}{subSqFine(j)}(k) = x(counter);
                                
                            end
                            
                            counter = counter + 1;
                            
                        end
                        
                    end
                    
                    %% Test -- Update total fix-ups and iterations
                    tot_fixups = tot_fixups + 1;
                    tot_iters = tot_iters + num_iters;
                
                end
                %}
                
            end
            
        end
        
    end
    
    %% Update coarse and fine information
    if orderFine < orderTo
        
        % Coarse quadrature order
        orderCoarse = orderCoarse + 1;
        
        % Coarse quadrature
        quadCoarse = GenQuad(orderCoarse);
        
        % Coarse quadrature solution
        % Note: Have to rearrange to match solFrom format
        numSubSqCoarse = size(quadCoarse, 1) / (3 * 4);
        counter = 1;
        for i = 1 : 3
            for j = 1 : numSubSqCoarse
                for k = 1 : 4
                    solCoarse(counter) = solFine{i}{j}(k);
                    counter = counter + 1;
                end
            end
        end
        
        % Fine quadrature order
        orderFine = orderFine + 1;
        
        % Fine quadrature
        quadFine = GenQuad(orderFine);
        
    else
        
        % Fine quadrature order
        orderFine = orderFine + 1;
        
    end
    
end

fprintf('Total number of fix-ups: %i \n', tot_fixups);
fprintf('Average number of iterations per fix-up: %e \n', tot_iters / tot_fixups);

end