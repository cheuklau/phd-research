function [iters, solution] = ft(c, A, b, eps1, eps2, eps3, bfs)
%
%	Solves:  minimize cx	subject to Ax <= b & x >= 0

%	m		number of rows in A
%	n		number of columns in A
%	B_indices	vector of columns in A comprising the solution basis
%	V_indices	vector of columns in A not in solution basis

[m n] = size(A);

B_indices = find(bfs);

V_indices = find(ones(1,n) - abs(sign(bfs)));

ft_nnz = zeros(5000,2);

refactors=0;

iters=0;

while 1 == 1
    
    refactors = refactors + 1;
    
    %	L		lower triangular factor of the basis
    %	U		upper triangular factor of the basis
    %	R		matrix of R vectors
    %	numR		number of entries in R
    %	Q		column permutation matrix
    
    B_indices = B_indices(:, colamd(A(:, B_indices)));
    
    [L U pt] = lu(A(:, B_indices));
    
    P = [1 : m] * pt';
    
    Q = [1 : m];
    
    R = zeros(m, m + 1);
    
    numR = 0;
    
    % Simplex method loops continuously until solution is found or discovered
    % to be impossible.
    
    ok = 1;
    
    while (nnz(R) < (.65 * m) & ok)
        
        iters = iters + 1;
        
        ft_nnz(iters, 1) = nnz(A(:, B_indices));
        ft_nnz(iters, 2) = nnz(L) + nnz(U) + nnz(R) - numR;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 1
        %	compute B^-1
        %	We won't explicitly calculate it here, but it will be solved when
        %	necessary to take advantage of matrix structures:
        %	B    = L * R(1) * ... * R(t) * U * Q'
        %	B^-1 = Q * U^-1 * Q * R(t)^-1 * ... * R(1)^-1 * L^-1
        
        %	Qinv		reverse column permutation of Q
        
        Qinv(Q) = [1 : m];
        
        Pinv(P) = [1 : m];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 2
        %	compute d = B^-1 * b = B \ b
        
        %	dtemp		accumulating value of d during R computation
        %	z		loop variable
        %	d		current solution vector
        
        dtemp = L \ b(P);
        
        for z = 1 : numR
            
            rtemp = [-R(z, 2 : R(z, 1)) 1 -R(z, R(z, 1) + 2 : m + 1)];
            
            dtemp(R(z, 1)) = rtemp * dtemp;
            
        end
        
        dtemp = U(Q, :) \ dtemp(Q, :);
        
        d = dtemp(Qinv, :);
        
        if (norm(A(:, B_indices) * d - b, inf) > eps3)
            
            ok = 0;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 3/Step 4/Step 5
        %	compute c_tilde = c_V - c_B * B^-1 * V
        
        %	ctemp		accumulating value during R computation
        %	c_tilde		modified cost vector
        
        c_tilde = zeros(1, n);
        
        ctemp = L \ A(P, V_indices);
        
        for z = 1 : numR
            
            rtemp = [-R(z, 2 : R(z, 1)) 1 -R(z, R(z, 1) + 2 : m + 1)];
            
            ctemp(R(z, 1), :) = rtemp * ctemp;
        
        end
        
        ctemp = U(Q, :) \ ctemp(Q, :);
        
        c_tilde(:, V_indices) = c(:, V_indices) - c(:, B_indices) * ctemp(Qinv, :);
        
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
            
            solution = zeros(n, 1);
            
            solution(B_indices, :) = d;
            
            return;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 8
        %	otherwise, compute w = B^-1 * a[j]
        
        %	wtemp		accumulating value of w during R computation
        %	w		relative weight (vector) of column entering the basis
        
        wtemp = L \ A(P, j);
        
        for z = 1 : numR
            
            rtemp = [-R(z, 2 : R(z, 1)) 1 -R(z, R(z, 1) + 2 : m + 1)];
            
            wtemp(R(z, 1)) = rtemp * wtemp;
            
        end
        
        wtemp = U(Q, :) \ wtemp(Q, :);
        
        w = wtemp(Qinv, :);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 9
        %	compute i s.t. w[i]>0 and d[i]/w[i] is minimized
        
        %	mn		minimum of d[i]/w[i] when w[i] > 0
        %	i		row corresponding to mn -- detemines outgoing column
        %	k		temporary storage variable
        
        zz = find (w > eps1)';
        
        [~, ii] = min (d(zz) ./ w (zz));
        
        i = zz(ii(1));
        
        k = B_indices(i);
        
        B_indices(i) = j;
        
        V_indices(j == V_indices) = k;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Update U %%%%%%%	** Forrest-Tomlin **
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        %% U is a row-permuted upper-triangular matrix
        %% Q*U is the proper upper-triangular matrix
        
        %	uc		spiked column of U resulting from new basis column
        %	ur		row corresponding to uc
        
        uc = Qinv(i);
        
        ur = Q(uc);
        
        %	r		new entry in R that eliminates the nonzero values
        %			in ur
        
        r = [zeros(1, uc) U(ur, uc + 1 : m)];
        
        r = r / U;
        
        numR = numR + 1;
        
        R(numR, :) = [ur r];
        
        % We already know the effects of r, so there is no need to explicitly
        % calculate them.
        
        if (uc < m)
            
            U(ur, uc + 1 : m) = zeros(1, m - uc);
            
        end
        
        % Update the incoming column according to the operations already done on U
        
        incoming = L \ A(P, j);
        
        for z = 1 : numR
            
            rtemp = [-R(z, 2 : R(z, 1)) 1 -R(z, R(z, 1) + 2 : m + 1)];
            incoming(R(z, 1)) = rtemp * incoming;
            
        end
        
        if (abs(incoming(ur) < eps2))
            
            ok = 0;
            
        end
        
        U = [U(:, [1 : uc - 1 uc + 1 : m]) incoming];
        
        % And then update Q
        Q = Q(:, [1 : uc - 1 uc + 1 : m uc]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %	Step 10
        %	REPEAT
        
    end
    
end