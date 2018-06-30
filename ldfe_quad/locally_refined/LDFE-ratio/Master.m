% Master file for cube-based LDFE quadrature sets
close force all;
clear;
clc;

% Refinement level
initialRef = 1; % Initial uniform refinement
localRef   = 3; % Number of local refinements

% Write to data file
fid = fopen('LDFEQUADRATIOLOCAL0.DAT', 'at');

% Go through each local refinement level
for iRef = initialRef + 1 : localRef
    
    fprintf('Currently on local refinement level: %i \n', iRef);
    
    % Square property
    squareInfo = SquareProp(initialRef, iRef, 0.0);
    
    % LDFE-center
    quadrature = LDFEratio(squareInfo, 4);    
    
    % Print quadrature for PDT
    PrintPDT(squareInfo, quadrature, fid);
   
end

% Close PDT quadrature file
fclose(fid);