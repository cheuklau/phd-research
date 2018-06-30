% Preliminary
close force all;
clear;
clc;

% Refinement level
order = 5;

% Initialize relative error
intError = [];

% Write to data file
fid = fopen('QUADLDFEQUAD.dat', 'at');

% Go through each refinement level
counter = 1;

for iRef = 0 : 1 : order
    
    % Indirect refinement order
    orderUse = 2 ^ iRef - 1;    
        
    fprintf('Currently on refinement level: %i \n', orderUse);
            
    % Square property
    squareInfo = SquareProp(orderUse);
    
    % LDFE-ratio for first guess ratios
    [~, rhoInitAll] = genFirstRatios(squareInfo, 2);

    % LDFE-sa
    quadrature = LDFEsa(squareInfo, rhoInitAll);
    
    % Print quadrature for PDT
    PrintPDT(quadrature, orderUse, fid);
      
    % Update plot counter
    counter = counter + 1;
   
end

% Close PDT quadrature file
fclose(fid);