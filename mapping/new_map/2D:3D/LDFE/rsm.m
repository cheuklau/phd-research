% c    = cost vector of the linear program
% A    = coefficient matrix of the linear program
% b    = constant right-hand side vector of the linear program
% bfs  = initial basic feasible solution
% eps1 = zero tolerances

function [iters, solution] = rsm (c, A, b, eps1, bfs)
%
%	Solves:  minimize cx	subject to Ax <= b & x >= 0
%
%   From "Introduction to Linear Programming, Applications and Extenstions"
%      by Richard B. Darst
%      page 101
%
%	m		number of rows in A
%	n		number of columns in A
%	B_indices	vector of columns in A comprising the solution basis
%	V_indices	vector of columns in A not in solution basis

[~, n]    = size(A);

B_indices = find(bfs);

V_indices = find(ones(1, n) - abs(sign(bfs)));

rsm_nnz   = zeros(5000, 2);

% Simplex method loops continuously until solution is found or discovered
% to be impossible.

iters = 0;

while 1 == 1
    
    iters = iters + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Step 1
    %	compute B^-1
    %
    %	Binv		inverse of the basis (directly computed)
    
    if size(A(:,B_indices), 1) ~= size(A(:, B_indices), 2)
        
        disp('here');
        
    end
    
    Binv = inv(A(:, B_indices));
    
    rsm_nnz(iters, 1) = nnz(A(:, B_indices));
    
    rsm_nnz(iters, 2) = nnz(Binv);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Step 2
    %	compute d = B^-1 * b
    %
    %	d		current solution vector
    
    d = Binv * b;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Step 3/Step 4/Step 5
    %	compute c_tilde = c_V - c_B * B^-1 * V
    %
    %	c_tilde		modified cost vector
    
    c_tilde = zeros(1, n);
    
    c_tilde(:, V_indices) = c(:, V_indices) - c(:, B_indices) * Binv * A(:, V_indices);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Step 6
    %	compute j s.t. c_tilde[j] <= c_tilde[k] for all k in V_indices
    %
    %   cj minimum cost value (negative) of non-basic columns
    %   j column in A corresponding to minimum cost value
    
    [cj j]=min(c_tilde);
    
    %disp(cj);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 7
    % if cj >= 0 , then we're done -- return solution which is optimal
    
    if cj >= -eps1
        
        solution = zeros(n,1);
        
        solution(B_indices,:) = d;
                
        return;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Step 8
    %	compute w = B^-1 * a[j]
    %
    %	w		relative weight (vector) of column entering the basis
    
    w = Binv * A(:,j);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Step 9
    %	compute i s.t. w[i]>0 and d[i]/w[i] is a smallest positive ratio
    %	swap column j into basis and swap column i out of basis
    %
    %	mn		minimum of d[i]/w[i] when w[i] > 0
    %	i		row corresponding to mn -- determines outgoing column
    %	k		temporary storage variable
        
    zz = find (w > eps1)';
    
    [~, ii] = min (d(zz) ./ w (zz)) ;
    
    i = zz(ii(1)) ;
    
    if (i == 0)
        
        error ('System is unbounded.');
        
    end
    
    k = B_indices(i);
    
    B_indices(i) = j;
    
    V_indices(j == V_indices) = k;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Step 10
    %	REPEAT
    
end