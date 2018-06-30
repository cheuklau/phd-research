function PlotSol(quadFrom, solFrom, quadTo, solTo, functionNumber)

% Preliminary
figure

% Step function with discontinuity within an LDFE range
numDiv = 1000;
if functionNumber == 1
    
    meshsize_1 = linspace(-1, -1 / 3, numDiv);
    y_1 = ones(1, 1000);
    
    meshsize_2 = linspace(-1 / 3, 1, numDiv);
    y_2 = zeros(1, numDiv);
    
    meshsize = [meshsize_1, meshsize_2];
    y = [y_1, y_2];
    
    % Step function with discontinuity at mu=0
elseif functionNumber == 2
    
    meshsize_1 = linspace(-1, 0, numDiv);
    y_1 = ones(1, numDiv);
    
    meshsize_2 = linspace(0, 1, numDiv);
    y_2 = zeros(1, numDiv);
    
    meshsize = [meshsize_1, meshsize_2];
    y = [y_1, y_2];
    
    % Smooth linear function
elseif functionNumber == 3
    
    meshsize = linspace(-1, 1, numDiv);
    y = zeros(numDiv, 1);
    for i = 1 : numDiv
        y(i) = 1 + meshsize(i);
    end
    
    % Exponential function with discontinuity at mu=0
elseif functionNumber == 4
    
    meshsize = linspace(-1, 1, numDiv);
    y = zeros(numDiv, 1);
    for i = 1 : numDiv
        y(i) = exp(-1 / meshsize(i));
    end
    
    % Preliminary exam problem
elseif functionNumber == 5
    
    % Problem parameters
    S       = 1.0; % Source strength
    sigma_t = 0.1; % Total cross section
    delta_x = 0.1; % Source thickness
    
    meshsize = linspace(-1, 1, numDiv);
    y = zeros(numDiv, 1);
    for i = 1 : numDiv
        y(i) = (S / sigma_t) * ...
            (1 - exp(-1 * sigma_t * delta_x / (2 * abs(meshsize(i)))));
    end    
end
analytical = plot(meshsize, y, 'Color', 'k');
set(analytical, 'LineWidth', 1.25);
hold on

% Quadrature mapping from information
muFrom         = quadFrom.mu;
numRegionsFrom = size(muFrom, 2);
solFromPlot    = zeros(numRegionsFrom * 2, 1);
muFromPlot     = zeros(numRegionsFrom * 2, 1);
counter = 1;
for i = 1 : numRegionsFrom
    for j = 1 : 2
        solFromPlot(counter) = solFrom{i}(j);
        muFromPlot(counter)  = muFrom{i}(j);
        counter = counter + 1;
    end
end

% Plot analytic solution
fromPlot = plot(muFromPlot, solFromPlot, 'k+');
set(fromPlot, 'MarkerSize', 8);
hold on

% Quadrature mapping to information
muTo         = quadTo.mu;
numRegionsTo = size(muTo, 2);
solToPlot    = zeros(numRegionsTo * 2, 1);
muToPlot     = zeros(numRegionsTo * 2, 1);
counter = 1;
for i = 1 : numRegionsTo
    for j = 1 : 2
        solToPlot(counter) = solTo{i}(j);
        muToPlot(counter)  = muTo{i}(j);
        counter = counter + 1;
    end
end

% Plot mapped solution
toPlot = plot(muToPlot, solToPlot, 'ro');
set(toPlot, 'MarkerSize', 8, 'Color', 'r');
hold on

% Plot options
grid on
set(gca,'FontSize',14)
xlabel('\mu','FontSize',18)
ylabel('f(\mu)','FontSize',18)
if numRegionsFrom > numRegionsTo
   title(sprintf('Fine-to-coarse mapping from %i to %i directions with fix-up', numRegionsFrom * 2, numRegionsTo * 2), 'FontSize', 18)
else
   title(sprintf('Coarse-to-fine mapping from %i to %i directions with fix-up', numRegionsFrom * 2, numRegionsTo * 2), 'FontSize', 18)
end
legend('Analytical Solution', 'Original Solution', 'Mapped Solution', 'Location','Southeast')
 

end