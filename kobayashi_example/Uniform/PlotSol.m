function PlotSol(prob)

% Retrieve prob information
R = prob.R;
Q = prob.Q;
S = prob.S;
P = prob.P;

% Plot analytic solution
gridSize = 100;
gamma = linspace(0, pi / 2, gridSize);
theta = linspace(0, pi / 2, gridSize);
omegaXref = zeros(gridSize);
omegaYref = zeros(gridSize);
omegaZref = zeros(gridSize);
solRef = zeros(gridSize);
for i = 1 : gridSize
    for j = 1 : gridSize
        
        % Calculate directional cosines        
        omegaXref(i, j) = cos(gamma(i)) * sin(theta(j));
        omegaYref(i, j) = sin(gamma(i)) * sin(theta(j));
        omegaZref(i, j) = cos(theta(j));       
        
        % Distance from detector to origin
        d = sqrt(P(1) ^ 2 + P(2) ^ 2 + P(3) ^ 2);
        
        % Maximum angle between detector and source
        t_max = atan(R / d);
        
        % Angle between quadrature direction and source
        dot_prod = ...
            P(1) * omegaXref(i, j) + ...
            P(2) * omegaYref(i, j) + ...
            P(3) * omegaZref(i, j);
        
        quad_mag = ...
            sqrt(omegaXref(i, j) ^ 2 + ...
            omegaYref(i, j) ^ 2 + ...
            omegaZref(i, j) ^ 2);
        
        t_i = acos( dot_prod / (d * quad_mag));
        
        % Non-zero angular flux if quadrature direction passes through sphere
        if t_i < t_max
            
            % Calculate chord length
            alpha = pi / 2 - t_i;
            
            r = d / tan(alpha);
            
            beta = asin(sin(alpha) * (r / R));
            
            delta = pi - beta - alpha;
            
            chord_1 = sin(delta) * r / sin(beta);
            
            gamma_2 = pi  - alpha;
            
            xi = asin(r * sin(gamma_2) / R);
            
            chi = pi - gamma_2 - xi;
            
            chord_2 = sin(chi) * r / sin(xi);
            
            chord = chord_1 + chord_2;
            
            % Calculate angular flux solution
            solRef(i, j) = (Q / S) * (1 - exp(-1 * S * chord));
            
            % Zero angular flux if quadrature flux does not pass through sphere
        else
            
            solRef(i, j) = 0;
            
        end
        
        % Adjust magnitude of plotted solution
        omegaXref(i, j) = solRef(i, j) * omegaXref(i, j);
        omegaYref(i, j) = solRef(i, j) * omegaYref(i, j);
        omegaZref(i, j) = solRef(i, j) * omegaZref(i, j);
        
    end
end
mesh(omegaXref, omegaYref, omegaZref, 'facecolor', 'none', 'EdgeColor', 'k');
hold on

title('Analytic solution for spherical source prob','FontSize', 18);
xlabel('\mu', 'FontSize', 18);
ylabel('\eta', 'FontSize', 18);
zlabel('\xi', 'FontSize', 18);

end