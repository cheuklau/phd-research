% Master file for cube-based LDFE quadrature sets
close force all;
clear;
clc;

% Refinement level
initialRef = 1; % Initial uniform refinement
localRef   = 3; % Number of local refinements

% Write to data file
fid = fopen('TEST.DAT', 'at');

% Go through each local refinement level
for iRef = initialRef + 1 : localRef
    
    fprintf('Currently on local refinement level: %i \n', iRef);
    
    % Square property
    squareInfo = SquareProp(initialRef, iRef);
    
    % Initial ratios
    [rhoInitAll1, rhoInitAll2, rhoInitAll3] = initialRatios(squareInfo, 2);
    
    % LDFE-center
    quadrature = LDFEsa(squareInfo, rhoInitAll1, rhoInitAll2, rhoInitAll3);    
    
    % Print quadrature for PDT
    PrintPDT(squareInfo, quadrature, fid);
   
end

% Close PDT quadrature file
fclose(fid);