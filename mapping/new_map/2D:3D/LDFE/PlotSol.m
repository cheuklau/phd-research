function PlotSol(quadTo, solTo, functionNumber, orderFrom, orderTo, problem, limits, useFixUp, factor_per)

%% Plot analytic solution

% Create new figure
figure

% Set grid size
if functionNumber == 3
    
    % Grid size
    gridSize = 50;
    
    % Retrieve problem information
    R = problem.source_radius;
    Q = problem.source_strength;
    S = problem.source_sigma_t;
    P = problem.detector_pos;
    
    % Distance from detector to origin
    d = sqrt(P(1) ^ 2 + P(2) ^ 2 + P(3) ^ 2);
    
    % Maximum angle between detector and source
    t_max = atan(R / d);
    
elseif functionNumber == 4 || functionNumber == 5
    
    % Grid size
    gridSize = 50;
    
    % Gamma and theta limits
    min_gamma = problem.min_gamma;
    max_gamma = problem.max_gamma;
    min_theta = problem.min_theta;
    max_theta = problem.max_theta;
    
else
    
    % Grid size
    gridSize = 25;
    
end

% Gamma and theta discretization
gamma = linspace(0, pi / 2, gridSize);
theta = linspace(0, pi / 2, gridSize);

% Initialize storage
solRef    = zeros(gridSize);
gammaPlot = zeros(gridSize);
thetaPlot = zeros(gridSize);

% Go through each discretization in gamma
for i = 1 : gridSize
    
    % Go through each discretization in theta
    for j = 1 : gridSize
        
        % Store angles
        gammaPlot(i, j) = gamma(i);
        thetaPlot(i, j) = theta(j);
        
        % Calculate directional cosines
        omegaX = cos(gamma(i)) * sin(theta(j));
        omegaY = sin(gamma(i)) * sin(theta(j));
        omegaZ = cos(theta(j));
        
        % Linear function
        if functionNumber == 1
            
            solRef(i, j) = 1 + omegaX + omegaY + omegaZ;
            
            % 5th-order function
        elseif functionNumber == 2
            
            solRef(i, j) = 1 + omegaX ^ 5 + omegaY ^ 5 + omegaZ ^ 5;
            
            % Spherical source function
        elseif functionNumber == 3
            
            % Angle between quadrature direction and source
            dot_prod = P(1) * omegaX + P(2) * omegaY + P(3) * omegaZ;
            
            t_i = acos(dot_prod / d);
            
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
            
            % Delta-gamma, delta-theta patch function
        elseif functionNumber == 4
            
            % Calculate theta and gamma of current direction
            thetaTmp = acos(omegaZ);
            gammaTmp = acos(omegaX / sin(thetaTmp));
            
            % Determine if current direction is within patch of interest
            if thetaTmp > min_theta - factor_per && thetaTmp < max_theta + factor_per && ...
                    gammaTmp > min_gamma - factor_per && gammaTmp < max_gamma + factor_per
                
                solRef(i, j) = 1;
                
            else
                
                solRef(i, j) = 0;
                
            end
            
            
        elseif functionNumber == 5
            
            % Calculate theta and gamma of current direction
            thetaTmp = acos(omegaZ);
            gammaTmp = acos(omegaX / sin(thetaTmp));
            
            % Determine if current direction is within patch of interest
            if thetaTmp > pi / 8 && thetaTmp < pi / 4 && ...
                    gammaTmp > pi / 8 && gammaTmp < 3 * pi / 8
                
                solRef(i, j) = 1 + factor_per;
                
            elseif thetaTmp > pi / 4 && thetaTmp < 3 * pi / 8 && ...
                    gammaTmp > pi / 8 && gammaTmp < 3 * pi / 8
                
                solRef(i, j) = 1 - factor_per;
                
            else
                
                solRef(i, j) = 0;
                
            end
            
        end
        
    end
    
end

% Analytic solution plot
surf(gammaPlot, thetaPlot, solRef, 'facecolor', 'none', 'EdgeColor', 'k');
hold on
%{
%title('Analytic solution for spherical source problem','FontSize', 18);
xlabel('\gamma', 'FontSize', 18);
ylabel('\theta', 'FontSize', 18);
zlabel('f(\Omega)', 'FontSize', 18);
%}
%% Analytic solution for quadrature you are mapping to

% Initialize storage
numDirsTo = size(quadTo, 1);
solRefTo  = zeros(numDirsTo, 1);
gammaPlot = zeros(numDirsTo, 1);
thetaPlot = zeros(numDirsTo, 1);

% Go through each direction you are mapping to
for i = 1 : numDirsTo
    
    % Store gamma and theta
    thetaPlot(i) = acos(quadTo{i}(3));
    gammaPlot(i) = acos(quadTo{i}(1) / sin(thetaPlot(i)));
    
    % Linear function
    if functionNumber == 1
        
        solRefTo(i) = 1 + quadTo{i}(1) + quadTo{i}(2) + quadTo{i}(3);
        
        % 5th-order function
    elseif functionNumber == 2
        
        solRefTo(i) = 1 + quadTo{i}(1) ^ 5 + quadTo{i}(2) ^ 5 + quadTo{i}(3) ^ 5;
        
        % Spherical source function
    elseif functionNumber == 3
        
        % Angle between quadrature direction and source
        dot_prod = P(1) * quadTo{i}(1) + P(2) * quadTo{i}(2) + P(3) * quadTo{i}(3);
        
        t_i = acos(dot_prod / d);
        
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
        
    elseif functionNumber == 4
        
        % Determine if current direction is within patch of interest
        if thetaPlot(i) > min_theta - factor_per && thetaPlot(i) < max_theta + factor_per && ...
                gammaPlot(i) > min_gamma - factor_per && gammaPlot(i) < max_gamma + factor_per
            
            solRefTo(i) = 1;
            
        else
            
            solRefTo(i) = 0;
            
        end
        
    elseif functionNumber == 5
        
        % Determine if current direction is within patch of interest
        if thetaPlot(i) > pi / 8 && thetaPlot(i) < pi / 4 && ...
                gammaPlot(i) > pi / 8 && gammaPlot(i) < 3 * pi / 8
            
            solRefTo(i, j) = 1 + factor_per;
            
        elseif thetaPlot(i) > pi / 4 && thetaPlot(i) < 3 * pi / 8 && ...
                gammaPlot(i) > pi / 8 && gammaPlot(i) < 3 * pi / 8
            
            solRefTo(i, j) = 1 - factor_per;
            
        else
            
            solRefTo(i, j) = 0;
            
        end
        
    end
    
end

%% Generate error plot

% Initialize storage
numSubSq  = numDirsTo / (4 * 3);
errorTo   = zeros(numDirsTo, 1);
solToPlot = zeros(numDirsTo, 1);
counter   = 1;
psiMax    = limits.psiMax;

% Go through each face
for i = 1 : 3
    
    % Go through each sub-square
    for j = 1 : numSubSq
        
        % Go through each direction
        for k = 1 : 4
            
            % Store mapped solution to plot
            solToPlot(counter) = solTo{i}{j}(k);
            
            % Store error
            if functionNumber == 3 || functionNumber == 4 || functionNumber == 5
                
                errorTo(counter) = (solTo{i}{j}(k) - solRefTo(counter)) / psiMax;
                
            else
                
                errorTo(counter) = (solTo{i}{j}(k) - solRefTo(counter)) / solRefTo(counter);
                
            end
            
            % Update counter
            counter = counter + 1;
            
        end
        
    end
    
end

%% Set plot options

fprintf('min error: %f \n', min(errorTo));
fprintf('max error: %f \n', max(errorTo));

% Error plot options
caxis([min(errorTo) max(errorTo)]);
caxis([-0.8855 0.9616]);
scatter3(gammaPlot, thetaPlot, solToPlot, 80, errorTo, 'filled');

axis([0 pi/2 0 pi/2 -1 2]);
set(gca, 'FontSize', 18);
errorColorBar = colorbar;
set(errorColorBar, 'fontsize', 18);
% Axes labeles
ylabel(errorColorBar, 'Absolute Error', 'FontSize', 18);
xlabel('\gamma', 'FontSize', 24);
ylabel('\theta', 'FontSize', 24);
zlabel('f(\Omega)', 'FontSize', 24);

%{
% Quadrature information
numSubSqFrom = (orderUseFrom + 1) ^ 2 * 3;
numDirsFrom  = 4 * numSubSqFrom;

% Plot title
if useFixUp == 0
    
    title(sprintf('LDFE-SQ mapping without fix-up: %i to %i directions', ...
        numDirsFrom, numDirsTo),'FontSize', 20);
    
else
    
    title(sprintf('LDFE-SQ mapping with fix-up: %i to %i directions', ...
        numDirsFrom, numDirsTo),'FontSize', 20);
    
end
%}

end