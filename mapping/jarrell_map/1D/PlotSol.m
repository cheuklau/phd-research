function PlotSol(quadFrom, solFrom, orderFrom, quadTo, solTo, ...
    orderTo, functionNumber)

% Preliminary
figure

% Step function with discontinuity within an LDFE range
if functionNumber == 1
    
    meshsize_1 = linspace(-1, -1 / 3, 1000);
    y_1 = ones(1, 1000);

    meshsize_2 = linspace(-1 / 3, 1, 1000);
    y_2 = zeros(1, 1000);

    meshsize = [meshsize_1, meshsize_2];
    y = [y_1, y_2];
    
    analytical = plot(meshsize, y, 'Color', 'k');
    set(analytical, 'LineWidth', 1.25);
    hold on 
    
    % Step function with discontinuity at mu=0
elseif functionNumber == 2
    
    meshsize_1 = linspace(-1, 0, 1000);
    y_1 = ones(1, 1000);

    meshsize_2 = linspace(0, 1, 1000);
    y_2 = zeros(1, 1000);

    meshsize = [meshsize_1, meshsize_2];
    y = [y_1, y_2];
    
    analytical = plot(meshsize, y, 'Color', 'k');
    set(analytical, 'LineWidth', 1.25);
    hold on 
    
    % Smooth linear function
elseif functionNumber == 3
    
    meshsize = linspace(-1, 1, 2000);
    y = zeros(2000, 1);
    for i = 1 : 2000
        y(i) = 1 + meshsize(i);
    end
    analytical = plot(meshsize, y, 'Color', 'k');
    set(analytical,'LineWidth',1.0);
    hold on
    
    % Exponential function with discontinuity at mu=0
elseif functionNumber == 4
    
    meshsize = linspace(-1, 1, 100);
    y = zeros(100, 1);
    for i = 1 : 100
        y(i) = exp(-1 / meshsize(i));
    end
    analytical = plot(meshsize, y, 'Color', 'k');
    set(analytical,'LineWidth',1.0);
    hold on

end

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
legend('Analytical Solution', ...
    sprintf('Original Solution - Order %i', orderFrom), ...
    sprintf('Mapped Solution - Order %i', orderTo), ...
     'Location','Southeast')
 
end