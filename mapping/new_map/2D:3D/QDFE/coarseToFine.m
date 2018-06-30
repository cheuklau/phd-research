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

%% TEST --- Total number of fix-ups and iterations
tot_iters  = 0;
tot_fixups = 0;

%% Go through each refinement level
while orderFine <= orderTo
    
    %% Sub-square properties
    
    % Number of fine sub-squares per face
    numSubSqFine    = size(quadFine, 1) / (3 * 9);
    numSubSqFineRow = numSubSqFine ^ 0.5;
    
    % Number of coarse sub-squares per face
    numSubSqCoarse = size(quadCoarse, 1) / (3 * 9);
    
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
            subSqFine(1) = deltaX + deltaY * numSubSqFineRow + 1;
            subSqFine(2) = subSqFine(1) + 1;
            subSqFine(3) = subSqFine(2) + 1;
            subSqFine(4) = subSqFine(1) + numSubSqFineRow;
            subSqFine(5) = subSqFine(4) + 1;
            subSqFine(6) = subSqFine(5) + 1;
            subSqFine(7) = subSqFine(4) + numSubSqFineRow;
            subSqFine(8) = subSqFine(7) + 1;
            subSqFine(9) = subSqFine(8) + 1;
            
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
                
                % Coarse direction index
                iDir = (iFace - 1) * numSubSqCoarse * 9 + (i - 1) * 9 + j;
                
                % Go through each corresponding fine sub-square
                for k = 1 : 9
                    
                    % Go through each fine direction
                    for m = 1 : 9
                        
                        % Current fine sub-square data line
                        lineUse = lineFine(k) + m - 1;
                        
                        % Directional cosines
                        mu = quadFine{lineUse}(1);
                        eta = quadFine{lineUse}(2);
                        xi = quadFine{lineUse}(3);
                        wt = quadFine{lineUse}(4);
                        
                        % Coarse basis evaluated at fine direction
                        basisValue = ...
                            quadCoarse{iDir}(5) + ...
                            quadCoarse{iDir}(6) * mu + ...
                            quadCoarse{iDir}(7) * eta + ...
                            quadCoarse{iDir}(8) * xi + ...
                            quadCoarse{iDir}(9) * (mu ^ 2 - eta ^ 2) + ...
                            quadCoarse{iDir}(10) * xi ^ 2 + ...
                            quadCoarse{iDir}(11) * mu * eta + ...
                            quadCoarse{iDir}(12) * mu * xi + ...
                            quadCoarse{iDir}(13) * eta * xi;
                        
                        % Add contribution to A-matrix
                        aMat(1, j) = aMat(1, j) + wt * basisValue;
                        aMat(2, j) = aMat(2, j) + wt * basisValue * mu;
                        aMat(3, j) = aMat(3, j) + wt * basisValue * eta;
                        aMat(4, j) = aMat(4, j) + wt * basisValue * xi;
                        aMat(5, j) = aMat(5, j) + wt * basisValue * (mu ^ 2 - eta ^ 2);
                        aMat(6, j) = aMat(6, j) + wt * basisValue * xi ^ 2;
                        aMat(7, j) = aMat(7, j) + wt * basisValue * mu * eta;
                        aMat(8, j) = aMat(8, j) + wt * basisValue * mu * xi;
                        aMat(9, j) = aMat(9, j) + wt * basisValue * eta * xi;
                        
                    end
                    
                end
                
            end
            
            % Initialize b-vector
            bvec = zeros(9, 1);
            
            % Go through each coarse direction
            for j = 1 : 9
                
                % Coarse direction index
                iDir = (iFace - 1) * numSubSqCoarse * 9 + (i - 1) * 9 + j;
                
                % Directional cosines
                mu = quadCoarse{iDir}(1);
                eta = quadCoarse{iDir}(2);
                xi = quadCoarse{iDir}(3);
                wt = quadCoarse{iDir}(4);
                sol = solCoarse(iDir);
                
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
            
            % Solve for psi-tilde values
            coarseSolTilde = aMat \ bvec;
            
            % Go over each fine sub-square
            for j = 1 : 9
                
                % Go over each fine direction
                for k = 1 : 9
                    
                    % Initialize values
                    solFine{iFace}{subSqFine(j)}(k) = 0;
                    
                    % Fine quadrature index
                    lineUse = lineFine(j) + k - 1;
                    
                    % Go over each coarse direction
                    for m = 1 : 9
                        
                        % Coarse quadrature index
                        iDir = (iFace - 1) * numSubSqCoarse * 9 + (i - 1) * 9 + m;
                        
                        % Directional cosines
                        mu = quadFine{lineUse}(1);
                        eta = quadFine{lineUse}(2);
                        xi = quadFine{lineUse}(3);
                        
                        % Mapped fine solution
                        solFine{iFace}{subSqFine(j)}(k) = ...
                            solFine{iFace}{subSqFine(j)}(k) + ...
                            coarseSolTilde(m) * ...
                            (quadCoarse{iDir}(5) + ...
                            quadCoarse{iDir}(6) * mu + ...
                            quadCoarse{iDir}(7) * eta + ...
                            quadCoarse{iDir}(8) * xi + ...
                            quadCoarse{iDir}(9) * (mu ^ 2 - eta ^ 2) + ...
                            quadCoarse{iDir}(10) * xi ^ 2 + ...
                            quadCoarse{iDir}(11) * mu * eta + ...
                            quadCoarse{iDir}(12) * mu * xi + ...
                            quadCoarse{iDir}(13) * eta * xi);
                        
                    end
                    
                end
                
            end
            
            % Update fine sub-square position
            deltaX = deltaX + 3;
            if deltaX == numSubSqFineRow
                deltaX = 0;
                deltaY = deltaY + 3;
            end
            
            %% Apply fix-up algorithm
            
            if useFixUp == 1
                
                %% Check mapped solution for extrema and calculate net change in scalar flux
                phi_net     = 0;
                phi_bad     = 0;
                mu_bad      = 0;
                eta_bad     = 0;
                xi_bad      = 0;
                mu2eta2_bad = 0;
                xi2_bad     = 0;
                mueta_bad   = 0;
                muxi_bad    = 0;
                etaxi_bad   = 0;
                numBadOrd = 0;
                good_ord  = zeros(1, 81);
                
                % Go through each fine sub-square
                for j = 1 : 9
                    
                    % Go through each direction of current sub-square
                    for k = 1 : 9
                        
                        lineUse = lineFine(j) + k - 1;
                        
                        mu  = quadFine{lineUse}(1);
                        eta = quadFine{lineUse}(2);
                        xi  = quadFine{lineUse}(3);
                        wt  = quadFine{lineUse}(4);
                        
                        if solFine{iFace}{subSqFine(j)}(k) > psiMax
                            
                            phi_net = phi_net + ...
                                wt * (solFine{iFace}{subSqFine(j)}(k) - psiMax);
                            
                            solFine{iFace}{subSqFine(j)}(k) = psiMax;
                            
                            phi_bad     = phi_bad     + wt * psiMax;
                            mu_bad      = mu_bad      + wt * mu  * psiMax;
                            eta_bad     = eta_bad     + wt * eta * psiMax;
                            xi_bad      = xi_bad      + wt * xi  * psiMax;
                            mu2eta2_bad = mu2eta2_bad + wt * (mu ^ 2 - eta ^ 2) * psiMax;
                            xi2_bad     = xi2_bad     + wt * xi ^ 2 * psiMax;
                            mueta_bad   = mueta_bad   + wt * mu  * eta * psiMax;
                            muxi_bad    = muxi_bad    + wt * mu  * xi * psiMax;
                            etaxi_bad   = etaxi_bad   + wt * eta * xi * psiMax;
                            
                            numBadOrd = numBadOrd + 1;
                            
                        elseif solFine{iFace}{subSqFine(j)}(k) < psiMin
                            
                            phi_net = phi_net + ...
                                wt * (solFine{iFace}{subSqFine(j)}(k) - psiMin);
                            
                            solFine{iFace}{subSqFine(j)}(k) = psiMin;
                            
                            phi_bad     = phi_bad     + wt * psiMin;
                            mu_bad      = mu_bad      + wt * mu * psiMin;
                            eta_bad     = eta_bad     + wt * eta * psiMin;
                            xi_bad      = xi_bad      + wt * xi * psiMin;
                            mu2eta2_bad = mu2eta2_bad + wt * (mu ^ 2 - eta ^ 2) * psiMin;
                            xi2_bad     = xi2_bad     + wt * xi ^ 2 * psiMin;
                            mueta_bad   = mueta_bad   + wt * mu * eta * psiMin;
                            muxi_bad    = muxi_bad    + wt * mu * xi * psiMin;
                            etaxi_bad   = etaxi_bad   + wt * eta * xi * psiMin;
                            
                            numBadOrd = numBadOrd + 1;
                            
                        else
                            
                            good_ord((j - 1) * 9 + k) = 1;
                            
                        end
                        
                    end
                    
                end
                
                %% Solve multi-objective problem
                if numBadOrd > 0
                    
                    %% Calculate 0th and 1st moments from coarse sub-square
                    
                    % Initialize storage
                    phi_total     = 0;
                    mu_total      = 0;
                    eta_total     = 0;
                    xi_total      = 0;
                    mu2eta2_total = 0;
                    xi2_total     = 0;
                    mueta_total   = 0;
                    muxi_total    = 0;
                    etaxi_total   = 0;

                    % Go through each coarse direction
                    for j = 1 : 9
                        
                        % Coarse quadrature index
                        iDir = (iFace - 1) * numSubSqCoarse * 9 + (i - 1) * 9 + j;
                        
                        mu  = quadCoarse{iDir}(1);
                        eta = quadCoarse{iDir}(2);
                        xi  = quadCoarse{iDir}(3);
                        wt  = quadCoarse{iDir}(4);
                        sol = solCoarse(iDir);
                        
                        % 0th moment for each coarse direction
                        phi_total = phi_total + wt * sol;
                        
                        % 1st moment for each coarse direction
                        mu_total  = mu_total  + wt * mu  * sol;
                        eta_total = eta_total + wt * eta * sol;
                        xi_total  = xi_total  + wt * xi  * sol;
                        
                        % 2nd moment for each coarse direction
                        mu2eta2_total = mu2eta2_total + wt * (mu ^ 2 - eta ^ 2) * sol;
                        xi2_total     = xi2_total     + wt * xi ^ 2 * sol;
                        mueta_total   = mueta_total   + wt * mu * eta * sol;
                        muxi_total    = muxi_total    + wt * mu * xi * sol;
                        etaxi_total   = etaxi_total   + wt * eta * xi * sol;
                        
                    end
                    
                    %% Determine good ordinate limits for positive net scalar flux
                    if phi_net > 0
                        
                        counter = 1;
                        
                        phi_temp = 0;
                        
                        psi_min_ind = zeros(1, 81);
                        
                        for j = 1 : 9
                            
                            for k = 1 : 9
                                
                                lineUse = lineFine(j) + k - 1;
                                
                                if good_ord(counter) == 1
                                    
                                    psi_min_ind(counter) = solFine{iFace}{subSqFine(j)}(k);
                                    
                                else
                                    
                                    psi_min_ind(counter) = -1;
                                    
                                end
                                
                                phi_temp = phi_temp + ...
                                    quadFine{lineUse}(4) * solFine{iFace}{subSqFine(j)}(k);
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                        counter = 1;
                        
                        psi_max_ind = zeros(1, 81);
                        
                        for j = 1 : 9 
                            
                            for k = 1 : 9
                                
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
                        
                        counter = 1;
                        
                        phi_temp = 0;
                        
                        psi_max_ind = zeros(1, 81);
                        
                        for j = 1 : 9
                            
                            for k = 1 : 9
                                
                                lineUse = lineFine(j) + k - 1;
                                
                                if good_ord(counter) == 1
                                    
                                    psi_max_ind(counter) = solFine{iFace}{subSqFine(j)}(k);
                                    
                                else
                                    
                                    psi_max_ind(counter) = -1;
                                    
                                end
                                
                                phi_temp = phi_temp + ...
                                    quadFine{lineUse}(4) * solFine{iFace}{subSqFine(j)}(k);
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                        counter = 1;
                        
                        psi_min_ind = zeros(1, 81);
                        
                        for j = 1 : 9 
                            
                            for k = 1 : 9
                                
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
                    numGoodOrd = 81 - numBadOrd;                                        
                    
                    %% A-matrix of constraint coefficients
                    
                    % Initialize A-matrix
                    Amat = zeros(9 + 2 * numGoodOrd, 25 + 4 * numGoodOrd);
                    
                    % 0th through 2nd moment objectives
                    counter = 1;
                    
                    for j = 1 : 9
                        
                        for k = 1 : 9
                            
                            lineUse = lineFine(j) + k - 1;
                            
                            if good_ord((j - 1) * 9 + k) == 1
                                
                                mu  = quadFine{lineUse}(1);
                                eta = quadFine{lineUse}(2);
                                xi  = quadFine{lineUse}(3);
                                wt  = quadFine{lineUse}(4);
                            
                                % 1st/2nd moment objectives
                                Amat(1, counter) = wt * mu * eta;
                                Amat(2, counter) = wt * mu * xi;
                                Amat(3, counter) = wt * eta * xi;
                                Amat(4, counter) = wt * xi ^ 2;
                                Amat(5, counter) = wt * (mu ^ 2 - eta ^ 2);
                                Amat(6, counter) = wt * mu;
                                Amat(7, counter) = wt * eta;
                                Amat(8, counter) = wt * xi;
                                Amat(9, counter) = wt;
                            
                                counter = counter + 1;
                            
                            end
                                                        
                        end
                        
                    end
                    
                    % Main diagonal for artificial and slack variables
                    for j = 1 : 9 + 2 * numGoodOrd
                        
                        Amat(j, 2 * numGoodOrd + 16 + j) = 1;
                        
                    end
                    
                    % Off diagonal for psi limits
                    for j = 1 : numGoodOrd
                        
                        Amat(9 + j, j) = 1;
                        Amat(9 + numGoodOrd + j, j) = 1;
                        
                    end
                    
                    % Off diagonals for 1st and 2nd moment deltas
                    for j = 1 : 8
                        
                        Amat(j, numGoodOrd + j) = -1;
                        Amat(j, numGoodOrd + 8 + j) = 1;
                        
                    end
                    
                    % Surplus variables for psi min limit
                    for j = 1 : numGoodOrd
                        
                        Amat(9 + j, numGoodOrd + 16 + j) = -1;
                        
                    end
                    
                    %% b-vector of constraint right-hand sides
                    bVec = zeros(9 + 2 * numGoodOrd, 1);
                    bVec(1) = mueta_total - mueta_bad;
                    bVec(2) = muxi_total - muxi_bad;
                    bVec(3) = etaxi_total - etaxi_bad;
                    bVec(4) = xi2_total - xi2_bad;
                    bVec(5) = mu2eta2_total - mu2eta2_bad;
                    bVec(6) = mu_total - mu_bad;
                    bVec(7) = eta_total - eta_bad;
                    bVec(8) = xi_total - xi_bad;
                    bVec(9) = phi_total - phi_bad;
                    
                    counter = 1;
                    
                    for j = 1 : 9
                        
                        for k = 1 : 9
                                                    
                            if good_ord((j - 1) * 9 + k) == 1
                            
                                bVec(9 + counter) = psi_min_ind((j - 1) * 9 + k);
                                
                                bVec(9 + numGoodOrd + counter) = psi_max_ind((j - 1) * 9 + k);
                                
                                counter = counter + 1;
                                
                            end
                            
                        end
                        
                    end                    
                    
                    %% c-vector of cost function coefficients
                    cVec = zeros(1, 25 + 4 * numGoodOrd);
                    for j = numGoodOrd + 1 : numGoodOrd + 16
                        cVec(j) = 1;
                    end
                    for j = 17 + 2 * numGoodOrd : 25 + 3 * numGoodOrd
                        cVec(j) = 100;
                    end
                    
                    %% Initial feasible solution for Simplex methods
                    if LPtype ~= 4
                        bfs = zeros(1, 25 + 4 * numGoodOrd);
                        bfs(2 * numGoodOrd + 17) = mueta_total   - mueta_bad;
                        bfs(2 * numGoodOrd + 18) = muxi_total    - muxi_bad;
                        bfs(2 * numGoodOrd + 19) = etaxi_total   - etaxi_bad;
                        bfs(2 * numGoodOrd + 20) = xi2_total     - xi2_bad;
                        bfs(2 * numGoodOrd + 21) = mu2eta2_total - mu2eta2_bad;
                        bfs(2 * numGoodOrd + 22) = mu_total - mu_bad;
                        bfs(2 * numGoodOrd + 23) = eta_total - eta_bad;
                        bfs(2 * numGoodOrd + 24) = xi_total - xi_bad;
                        bfs(2 * numGoodOrd + 25) = phi_total - phi_bad;  
                        
                        counter = 1;
                        
                        for j = 1 : 9
                            
                            for k = 1 : 9
                            
                                if good_ord((j - 1) * 9 + k) == 1
                                   
                                    bfs(2 * numGoodOrd + 25 + counter) = psi_min_ind((j - 1) * 9 + k);
                                    
                                    bfs(3 * numGoodOrd + 25 + counter) = psi_max_ind((j - 1) * 9 + k);
                                    
                                    counter = counter + 1;
                                    
                                end
                                
                            end
                            
                        end
                                                                       
                    end
                    
                    bVectest = bVec;
                    AMattest = Amat;
                    cVectest = cVec;
                    bfstest = bfs;
                    
                    %% Remove entries corresponding to zero right-hand side constraints
                    row_to_delete = [];
                    col_to_delete = [];
                    
                    for j = 1 : 9 + 2 * numGoodOrd
                        
                        if bVec(j) == 0
                            
                            row_to_delete = [row_to_delete j];
                            
                            % 1st/2nd moment deltas and artificial
                            if j < 9
                                
                                col_to_delete = [col_to_delete j+numGoodOrd j+numGoodOrd+8 j+2*numGoodOrd+16];
                                
                                % 0th moment artificial
                            elseif j < 10
                                
                                col_to_delete = [col_to_delete j+2*numGoodOrd+16];
                                
                                % Min limit surplus and artificial
                            elseif j < numGoodOrd + 10
                                
                                col_to_delete = [col_to_delete j+numGoodOrd+7 j+2*numGoodOrd+16];
                                
                                % Max limit slack
                            else
                                
                                col_to_delete = [col_to_delete j+2*numGoodOrd+16];
                                
                            end
                            
                        end
                        
                    end
                    
                    % Delete rows and columns
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
                    
                    %{
                    if min(x(1:numGoodOrd)) < 0
                        
                       pause;
                       
                    end
                    %}
                    
                    %% Store results into mapped solution
                    
                    % Initialize parameters
                    counter = 1;
                    
                    % Go through each fine sub-square
                    for j = 1 : 9
                        
                        % Go through each fine direction
                        for k = 1 : 9
                            
                            if good_ord((j - 1) * 9 + k) == 1
                                
                                % Set near zero solutions to zero
                                if abs(x(counter)) < 1e-10
                                
                                    solFine{iFace}{subSqFine(j)}(k) = 0;
                                
                                    % Store remaining solutions
                                else
                                
                                    solFine{iFace}{subSqFine(j)}(k) = x(counter);
                                
                                end
                            
                                % Increment counter
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
                phiInd     = zeros(9, 1);
                muInd      = zeros(9, 1);
                etaInd     = zeros(9, 1);
                xiInd      = zeros(9, 1);
                mu2eta2Ind = zeros(9, 1);
                xi2Ind     = zeros(9, 1);
                muetaInd   = zeros(9, 1);
                muxiInd    = zeros(9, 1);
                etaxiInd   = zeros(9, 1);
                
                % Go through each coarse direction
                for j = 1 : 9
                    
                    % Coarse quadrature index
                    iDir = (iFace - 1) * numSubSqCoarse * 9 + (i - 1) * 9 + j;
                    
                    mu  = quadCoarse{iDir}(1);
                    eta = quadCoarse{iDir}(2);
                    xi  = quadCoarse{iDir}(3);
                    wt  = quadCoarse{iDir}(4);
                    sol = solCoarse(iDir);
                    
                    % 0th moment for each coarse direction
                    phiInd(j) = wt * sol;
                    
                    % 1st moment for each coarse direction
                    muInd(j)  = wt * mu * sol;
                    etaInd(j) = wt * eta * sol;
                    xiInd(j)  = wt * xi * sol;
                    
                    % 2nd moment for each coarse direction
                    mu2eta2Ind(j) = wt * (mu ^ 2 - eta ^ 2) * sol;
                    xi2Ind(j)     = wt * xi ^ 2 * sol;
                    muetaInd(j)   = wt * mu * eta * sol;
                    muxiInd(j)    = wt * mu * xi * sol;
                    etaxiInd(j)   = wt * eta * xi * sol;
                    
                end
                
                % Total 1st moment
                totMu  = sum(muInd);
                totEta = sum(etaInd);
                totXi  = sum(xiInd);
                
                % Total 2nd moment
                totmu2eta2 = sum(mu2eta2Ind);
                totxi2     = sum(xi2Ind);
                totmueta   = sum(muetaInd);
                totmuxi    = sum(muxiInd);
                totetaxi   = sum(etaxiInd);
                
                %% Check mapped solutions for extrema
                
                % Initialize out-of-range ordinate counter
                numBadOrd  = 0;
                
                % Go over each fine sub-square
                for j = 1 : 9
                    
                    % Go over each fine direction
                    for k = 1 : 9
                        
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
                    Amat = zeros(179, 357);
                    
                    % 0th moment constraint and 1st/2nd moment objectives
                    counter = 1;
                    
                    for j = 1 : 9
                       
                        for k = 1 : 9
                           
                            lineUse = lineFine(j) + k - 1;
                            
                            mu  = quadFine{lineUse}(1);
                            eta = quadFine{lineUse}(2);
                            xi  = quadFine{lineUse}(3);
                            wt  = quadFine{lineUse}(4);
                            
                            % 1st/2nd moment objectives
                            Amat(1, counter) = wt * mu * eta;
                            Amat(2, counter) = wt * mu * xi;
                            Amat(3, counter) = wt * eta * xi;
                            Amat(4, counter) = wt * xi ^ 2;
                            Amat(5, counter) = wt * (mu ^ 2 - eta ^ 2);
                            Amat(6, counter) = wt * mu;
                            Amat(7, counter) = wt * eta;
                            Amat(8, counter) = wt * xi;
                            
                            % 0th moment constraints
                            Amat(8 + j, (j - 1) * 9 + k) = wt;
                            
                            counter = counter + 1;
                            
                        end
                        
                    end
                    
                    % Main diagonal for artificial and slack variables
                    for j = 1 : 179
                       
                        Amat(j, 178 + j) = 1;
                        
                    end
                    
                    % Off diagonal for psi limits
                    for j = 1 : 81
                       
                        Amat(17 + j, j) = 1;
                        Amat(98 + j, j) = 1;
                        
                    end
                    
                    % Off diagonals for 1st/2nd moment deltas
                    for j = 1 : 8
                       
                        Amat(j, 81 + j) = -1;
                        Amat(j, 89 + j) = 1;
                        
                    end
                    
                    % Surplus variables for psi min limit
                    for j = 1 : 81
                       
                        Amat(17 + j, 97 + j) = -1;
                        
                    end
                    
                    %% b-vector of constraint right-hand sides
                    bVec = zeros(179, 1);
                    bVec(1) = totmueta;
                    bVec(2) = totmuxi;
                    bVec(3) = totetaxi;
                    bVec(4) = totxi2;
                    bVec(5) = totmu2eta2;
                    bVec(6) = totMu;
                    bVec(7) = totEta;
                    bVec(8) = totXi;
                    for j = 1 : 9
                       bVec(8 + j) = phiInd(j);
                    end
                    for j = 1 : 81
                       bVec(17 + j) = psiMin;
                       bVec(98 + j) = psiMax;
                    end
                    
                    %% c-vector of cost function coefficients
                    cVec = zeros(1, 357);
                    for j = 82 : 97
                       cVec(j) = 1;
                    end
                    for j = 179 : 276
                       cVec(j) = 10;
                    end
                    
                    %% Initial feasible solution for Simplex methods
                    if LPtype ~= 4
                        bfs = zeros(1, 357);
                        bfs(179) = totmueta;
                        bfs(180) = totmuxi;
                        bfs(181) = totetaxi;
                        bfs(182) = totxi2;
                        bfs(183) = totmu2eta2;
                        bfs(184) = totMu;
                        bfs(185) = totEta;
                        bfs(186) = totXi;
                        for j = 1 : 9
                           bfs(186 + j) = phiInd(j);
                        end
                        for j = 1 : 81
                           bfs(195 + j) = psiMin;
                           bfs(276 + j) = psiMax;
                        end
                    end
                    
                    %% Remove entries corresponding to zero right-hand side constraints
                    row_to_delete = [];
                    col_to_delete = [];
                                       
                    for j = 1 : 179
                        
                        if bVec(j) == 0
                            
                            row_to_delete = [row_to_delete j];
                                                        
                            % 1st/2nd moment deltas and artificial
                            if j < 9
                                
                                col_to_delete = [col_to_delete 81+j 89+j 178+j];
                                
                                % 0th moment artificial
                            elseif j < 18
                                
                                tmp = j - 8;
                                col_to_delete = [col_to_delete 186+tmp];
                                
                                % Min limit surplus and artificial
                            elseif j < 99
                                
                                tmp = j - 17;
                                col_to_delete = [col_to_delete 97+tmp 195+tmp];
                                
                                % Max limit slack
                            else
                                
                                tmp = j - 98;
                                col_to_delete = [col_to_delete 276+tmp];
                                
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
                    
                    % Initialize parameters
                    counter = 1;
                    
                    % Go through each fine sub-square
                    for j = 1 : 9
                        
                        % Go through each fine direction
                        for k = 1 : 9
                            
                            % Set near zero solutions to zero
                            if abs(x(counter)) < 1e-10
                                
                                solFine{iFace}{subSqFine(j)}(k) = 0;
                                
                                % Store remaining solutions
                            else
                                
                                solFine{iFace}{subSqFine(j)}(k) = x(counter);
                                
                            end
                            
                            % Increment counter
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
        numSubSqCoarse = size(quadCoarse, 1) / (3 * 9);
        counter = 1;
        for i = 1 : 3
            for j = 1 : numSubSqCoarse
                for k = 1 : 9
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