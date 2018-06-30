function [iters, solution] = sbg (c, A, b, eps1, eps3, bfs)
%
%	Solves:  minimize cx	subject to Ax <= b & x>= 0

%	m		number of rows in A
%	n		number of columns in A
%	B_indices	vector of columns in A comprising the solution basis
%	V_indices	vector of columns in A not in solution basis

[m n] = size(A);
B_indices = find(bfs);
V_indices = find(ones(1,n) - abs(sign(bfs)));

sbg_nnz = zeros(5000,2);

% Simplex method loops continuously until solution is found or discovered
% to be impossible.

refactors=0;
iters=0;

while 1==1
    refactors = refactors + 1;
    
    %	L		lower triangular factor of the basis
    %	U		upper triangular factor of the basis
    %	pt		row permutation factor due to partial pivoting
    %			during LU factorization of the basis
    %	Q		column permutation vector
    
    B = A(:,B_indices);
    
    Qinv = colamd(B);
    [L U pt] = lu(B(:,Qinv));
    
    %	Qinv		reverse column permutation of Q
    
    Q(Qinv) = [1:m];
    
    %	Parray		array of row permutations that develop during
    %			the Bartels-Golub method
    
    eta_limit = 0.68*nnz(L);
    
    Parray = sparse(zeros(eta_limit, m));
    Parray(1,:) = ([1:m] * pt') - [1:m];
    numP = 1;
    
    %	eta_array	array of eta matrices that develop during the
    %			Bartels-Golub method
    
    numeta = 1;
    eta_array = zeros(eta_limit + 5, 3);
    
    eta_temp = eta_decomposition(L);
    eta_rows = size(eta_temp,1);
    eta_array(numeta+1:numeta+eta_rows, :) = eta_temp;
    numeta = numeta + eta_rows;
    
    ok = 1;
    
    while (numeta < eta_limit & ok)
        iters=iters+1;
        
        sbg_nnz(iters,1) = nnz(A(:,B_indices));
        sbg_nnz(iters,2) = nnz(L) + nnz(U) + numeta;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 1
        %	compute B^-1
        
        %	Rather, solve z=Bx when necessary to take advantage of U's
        %	structure.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 2
        %	compute d = B^-1 * b
        
        %	dtemp		accumulating value of d during eta matrix computation
        %	lp		counting variable to identify row permutations
        %	z		loop variable
        %	d		current solution vector
        
        dtemp = b;
        lp = 1;
        
        for z = 1:numeta
            if (eta_array(z,1) == 0)	% Time to permute
                dtemp = dtemp(Parray(lp,:) + [1:m], :);
                lp = lp + 1;
            else				% Still working on L^-1
                row = eta_array(z,2);
                col = eta_array(z,3);
                dtemp(row,:) = dtemp(row,:) + eta_array(z,1) * dtemp(col,:);
            end;
        end;
        
        d = U(:,Q) \ dtemp;
        
        if (norm(A(:,B_indices)*d - b, inf) > eps3)
            ok = 0;
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 3/Step 4/Step 5
        %	compute c_tilde = c_V - c_B * B^-1 * V
        
        %	ctemp		accumulating value during eta matrix computation
        %	c_tilde		modified cost vector
        
        ctemp = A(:,V_indices);
        lp = 1;
        
        for z = 1:numeta
            if (eta_array(z,1) == 0)	% Time to permute
                ctemp = ctemp(Parray(lp,:) + [1:m], :);
                lp = lp + 1;
            else				% Still working on L^-1
                row = eta_array(z,2);
                col = eta_array(z,3);
                ctemp(row,:) = ctemp(row,:) + eta_array(z,1) * ctemp(col,:);
            end;
        end;
        
        ctemp = U(:,Q) \ ctemp;
        
        c_tilde = zeros(1,n);
        c_tilde(:,V_indices) = c(:,V_indices) - c(:,B_indices) * ctemp;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 6
        %	compute j s.t. c_tilde[j] <= c_tilde[k] for all k in V_indices
        
        %	cj		minimum cost value (negative) of non-basic columns
        %	j		column in A corresponding to minimum cost value
        
        [cj j] = min(c_tilde);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 7
        %	if cj >= 0 , then we're done -- return solution which is optimal
        
        if cj >= -eps1
            solution = zeros(n,1);
            solution(B_indices,:) = d;
            return;
        end;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 8
        %	otherwise, compute w = B^-1 * a[j]
        
        %	wtemp		accumulating value of w during eta matrix computation
        %	w		relative weight (vector) of column entering the basis
        
        wtemp = A(:,j);
        lp = 1;
        
        for z = 1:numeta
            if (eta_array(z,1) == 0)	% Time to permute
                wtemp = wtemp(Parray(lp,:) + [1:m], :);
                lp = lp + 1;
            else				% Still working on L^-1
                row = eta_array(z,2);
                col = eta_array(z,3);
                wtemp(row,:) = wtemp(row,:) + eta_array(z,1) * wtemp(col,:);
            end;
        end;
        
        w = U(:,Q) \ wtemp;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 9
        %	compute i s.t. w[i]>0 and d[i]/w[i] is a smallest positive ratio
        
        %	mn		minimum of d[i]/w[i] when w[i] > 0
        %	i		row corresponding to mn -- detemines outgoing column
        %	k		temporary storage variable
        
        mn = inf;
        i = 0;
        
        zz = find (w > eps1)' ;
        [yy, ii] = min (d(zz) ./ w (zz)) ;
        i = zz(ii(1)) ;
        
        if (i == 0)
            error ('System is unbounded.');
        end;
        
        k = B_indices(i);
        B_indices(i) = j;
        V_indices(j == V_indices) = k;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Update U %%%%%%%	** Bartels-Golub **
        %%%%%%%%%%%%%%%%%%%%%%%%
        %
        %	P * U * Q = M(r) * M(r-1) * ... * M(1) * B
        %
        
        U = A(:,B_indices);
        U = U(:,Qinv);
        lp = 1;
        
        for z = 1:numeta
            if (eta_array(z,1) == 0)	% Time to permute
                U = U(Parray(lp,:) + [1:m],:);
                lp = lp + 1;
            else				% Still working on L^-1
                row = eta_array(z,2);
                col = eta_array(z,3);
                U(row,:) = U(row,:) + eta_array(z,1) * U(col,:);
            end;
        end;
        
        % U will be the same as the U that started this iteration with the exception
        % of a spike in the column that corresponds to the column that was replaced
        % in B (column i)
        
        % 	lside		left side and top of bump (column/row in U)
        %	bside		right side and bottom of bump (column/row in U)
        
        lside = Q(i);
        bside = max(find(U(:,lside)));
        
        %	Ptemp		determines entire row permutation for this
        %			iteration of Bartels-Golub
        
        Ptemp = [1:m];
        numeta = numeta + 1;
        if (numeta < eta_limit)
            eta_array = [eta_array; 0 0 0];
        end;  % if
        
        % (Reid's method rotates the columns and rows here)
        % reidrotate;
        % Q(Qinv) = [1:m];
        
        if (lside < bside)
            
            %% Now do the Hessenberg permutation that was proposed in the sparse Bartels-
            %% Golub algorithm
            
            U    = U   (:, [1:lside-1 lside+1:bside lside bside+1:m]);
            Qinv = Qinv(:, [1:lside-1 lside+1:bside lside bside+1:m]);
            Q(Qinv) = [1:m];
            
            % And factor the Hessenberg portion of U, updating eta_array and Parray as
            % necessary.
            
            %	newL		lower triangular factor of the bump
            %	newU		upper triangular factor of the bump
            %	newP		row permutation occurring during LU factorization
            %			of the bump
            
            [newL newU newP] = lu(U(lside:bside, lside:bside));
            
            % Fix U
            U(lside:bside, :) = newP * U(lside:bside, :);
            U(lside:bside,lside:bside) = newU;
            
            if (bside < m)
                U(lside:bside, bside+1:m) = newL \ U(lside:bside, bside+1:m);
            end;
            
            % Fix P (Parray)
            Ptemp(lside:bside) = Ptemp(lside:bside) * newP';
            
            % Fix L (eta_array)
            Ltemp = eye(m);
            Ltemp(lside:bside, lside:bside) = newL;
            new_etas = eta_decomposition(Ltemp);
            if (new_etas ~= [])
                eta_rows = size(new_etas,1);
                eta_array(numeta+1:numeta+eta_rows, :) = new_etas;
                numeta = numeta + eta_rows;
            end;	% if
        end;	% if
        
        Ptemp = Ptemp - [1:m];
        Parray(numP+1,:) = Ptemp;
        numP = numP + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 10
        %	REPEAT
        
    end;	% while
    
end;  %while