function PlotSol(quadTo, solTo, functionNumber, orderFrom, orderTo, problem)

% Plot analytic solution
if functionNumber == 3
    gridSize = 50;
else
    gridSize = 25;
end
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
        
        % Linear function
        if functionNumber == 1
            
            solRef(i, j) = 1 + omegaXref(i, j) + omegaYref(i, j) + omegaZref(i, j);
            
            % Quadratic function
        elseif functionNumber == 2
            
            solRef(i, j) = 1 + omegaXref(i, j) ^ 5 + omegaYref(i, j) ^ 5 + omegaZref(i, j) ^ 5;
            
            % Spherical source function
        elseif functionNumber == 3
            
            % Retrieve problem information
            R = problem.source_radius;
            Q = problem.source_strength;
            S = problem.source_sigma_t;
            P = problem.detector_pos;
            
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
            
            
        end
        
        % Adjust magnitude of plotted solution
        omegaXref(i, j) = solRef(i, j) * omegaXref(i, j);
        omegaYref(i, j) = solRef(i, j) * omegaYref(i, j);
        omegaZref(i, j) = solRef(i, j) * omegaZref(i, j);
        
    end
end
mesh(omegaXref, omegaYref, omegaZref, 'facecolor', 'none', 'EdgeColor', 'k');
hold on

% Analytic solution for quadrature you are mapping to
numDirsTo = size(quadTo, 1);
solRefTo = zeros(numDirsTo, 1);
for i = 1 : numDirsTo
    
    % Linear function
    if functionNumber == 1
        
        solRefTo(i) = 1 + quadTo{i}(1) + quadTo{i}(2) + quadTo{i}(3);
        
        % Quadratic function
    elseif functionNumber == 2
        
        solRefTo(i) = 1 + quadTo{i}(1) ^ 5 + quadTo{i}(2) ^ 5 + quadTo{i}(3) ^ 5;
        
    elseif functionNumber == 3
        
        % Angle between quadrature direction and source
        dot_prod = ...
            P(1) * quadTo{i}(1) + ...
            P(2) * quadTo{i}(2) + ...
            P(3) * quadTo{i}(3);
        
        quad_mag = ...
            sqrt(quadTo{i}(1) ^ 2 + ...
                 quadTo{i}(2) ^ 2 + ...
                 quadTo{i}(3) ^ 2);
             
        t_i = acos(dot_prod / (d * quad_mag));
        
        % Non-zero angular flux if quadrature direction passes through sphere
        if t_i < t_max
            
            % Calculate chord length
            alpha = pi / 2 - t_i;
            
            r = d / tan(alpha);
            
            beta = asin(sin(alpha) * (r / R));
            
            delta = pi - beta - alpha;
            
            chord_1 = sin(delta) * r / sin(beta);
            
            gamma = pi  - alpha;
            
            xi = asin(r * sin(gamma) / R);
            
            chi = pi - gamma - xi;
            
            chord_2 = sin(chi) * r / sin(xi);
            
            chord = chord_1 + chord_2;
            
            % Calculate angular flux solution
            solRefTo(i) = (Q / S) * (1 - exp(-1 * S * chord));
            
            % Zero angular flux if quadrature flux does not pass through sphere
        else
            
            solRefTo(i) = 0;
            
        end

    end
    
end

% Plot mapped solution with error pseudocolor
numSubSq = numDirsTo / (4 * 3);
omegaXto = zeros(numDirsTo, 1);
omegaYto = zeros(numDirsTo, 1);
omegaZto = zeros(numDirsTo, 1);
errorTo  = zeros(numDirsTo, 1);
counter = 1;
for i = 1 : 3
    for j = 1 : numSubSq
        for k = 1 : 4
            omegaXto(counter) = solTo{i}{j}(k) * quadTo{counter}(1);
            omegaYto(counter) = solTo{i}{j}(k) * quadTo{counter}(2);
            omegaZto(counter) = solTo{i}{j}(k) * quadTo{counter}(3);
            if solRefTo(counter) ~= 0
                errorTo(counter) = ...
                    (solTo{i}{j}(k) - solRefTo(counter)) / solRefTo(counter);
            else
                errorTo(counter) = 0;
            end
            counter = counter + 1;
        end
        
    end
end
caxis([min(errorTo) max(errorTo)]);
scatter3(omegaXto, omegaYto, omegaZto, 80, errorTo, 'filled');
errorColorBar = colorbar;
ylabel(errorColorBar, 'Relative Error', 'FontSize', 18);
xlabel('\mu', 'FontSize', 18);
ylabel('\eta', 'FontSize', 18);
zlabel('\xi', 'FontSize', 18);
orderUseFrom = 2 ^ orderFrom - 1;
numSubSqFrom = (orderUseFrom + 1) ^ 2 * 3;
numDirsFrom  = 4 * numSubSqFrom;
orderUseTo = 2 ^ orderTo - 1;
numSubSqTo = (orderUseTo + 1) ^ 2 * 3;
numDirsToTest = 4 * numSubSqTo;
if numDirsToTest ~= numDirsTo
    error('something wrong with your direction calculation');
end
title(sprintf('Jarrell/Adams Mapping 2: %i to %i directions (%i to %i SQ)', numDirsFrom, numDirsTo, numSubSqFrom, numSubSqTo),'FontSize', 18);

%{
if functionNumber == 1
    title(sprintf('Order %i to Order %i : f(\\Omega)=1+\\mu', orderFrom, orderTo),'FontSize', 18);
elseif functionNumber == 2
    title(sprintf('Order %i to Order %i : f(\\Omega)=1+\\mu^{2}', orderFrom, orderTo),'FontSize', 18);
elseif functionNumber == 3
    title(sprintf('Order %i to Order %i : Spherical Source Problem', orderFrom, orderTo),'FontSize', 18);
end
%}

end