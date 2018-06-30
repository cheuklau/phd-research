function analytic_sol = analytic_calc(v_top, v_bot, u_tip_tilde, flip, h, b, fcn)

% Initialize solution
analytic_sol = 0;

% Gauss order
N = 128;

% Divide v using N-point Gauss quadrature
[fv, wv] = lgwt(N, v_bot, v_top);

% Go over each v quadrature point
for i = 1 : N
   
    % Determine u bounds for rightside-up triangle
    if flip == 0
       
        u_min = (fv(i) - (v_top - (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));
       
        u_max = (-1 * fv(i) + (v_top + (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));
        
        % Determine u bounds for upside-down triangle
    else
        
        u_min = (-1 * fv(i) + (v_bot + (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));
        
        u_max = (fv(i) - (v_bot - (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));
     
    end
    
    % Divide u using N-point Gauss quadrature
    [fu, wu] = lgwt(N, u_min, u_max);
    
    % Go through each u quadrature point
    for j = 1 : N
        
        % Jacobian
        jac = h * (h - fv(i)) / (2 * (((fu(j) ^ 2 + 1) / 2) * (h - fv(i)) ^ 2 + fv(i) ^ 2) ^ (3 / 2));
        
        % omega-x * omega-y
        if fcn == 1
            
            r = sqrt((fu(j) ^ 2 + 1) * (1 / 2) * (h - fv(i)) ^ 2 + fv(i) ^ 2);
            omegaX = ((1 - fu(j)) / 2) * (h - fv(i)) / r;
            omegaY = ((1 + fu(j)) / 2) * (h - fv(i)) / r;
            analytic_sol = analytic_sol + jac * omegaX * omegaY * wv(i) * wu(j); 
            
            % omega-x ^ 3 * omega-y * omega-z
        elseif fcn == 2
            
            r = sqrt((fu(j) ^ 2 + 1) * (1 / 2) * (h - fv(i)) ^ 2 + fv(i) ^ 2);
            omegaX = ((1 - fu(j)) / 2) * (h - fv(i)) / r;
            omegaY = ((1 + fu(j)) / 2) * (h - fv(i)) / r;
            omegaZ = fv(i) / r;
            analytic_sol = analytic_sol + jac * omegaX ^ 3 * omegaY * omegaZ * wv(i) * wu(j); 
            
            % omega-x ^ 3 * omgea-y ^ 6 * omega-z ^ 15
        elseif fcn == 3
            
            r = sqrt((fu(j) ^ 2 + 1) * (1 / 2) * (h - fv(i)) ^ 2 + fv(i) ^ 2);
            omegaX = ((1 - fu(j)) / 2) * (h - fv(i)) / r;
            omegaY = ((1 + fu(j)) / 2) * (h - fv(i)) / r;
            omegaZ = fv(i) / r;
            analytic_sol = analytic_sol + jac * omegaX ^ 3 * omegaY ^ 6 * omegaZ ^ 15 * wv(i) * wu(j);
            
        end                        
        
    end

end

end