% Master file for cube-based QDFE quadrature sets
close force all;
clear;
clc;

% Refinement level 
initialRef = 0; % Initial uniform refinement
localRef   = 3; % Number of local refinement

% Write to data file
fid = fopen('QDFEQUADSALOCAL60FEATHER.DAT', 'at');

for iRef = initialRef + 1 : localRef
    
    fprintf('Currently on refinement level: %i \n', iRef);

    % Square property
    squareInfo = SquareProp(initialRef, iRef, 0.6);
    
    % LDFE-center
    quadrature = QDFEsa(squareInfo);    

    % Print quadrature for PDT
    PrintPDT(squareInfo, quadrature, fid);
    
end

% Close PDT quadrature file
fclose(fid);