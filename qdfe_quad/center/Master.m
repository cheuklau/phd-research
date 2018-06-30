% Master file for cube-based QDFE quadrature sets

close force all;
clear;
clc;

% Refinement level 
nRef = 10;

% Initialize relative error
intError = [];

% Write to data file
fid = fopen('QDFEcenter0to10.DAT', 'at');

% Go through each refinement level
counter = 1;

for iRef = 0 : 1 : nRef
    
    fprintf('Currently on refinement level: %i \n', iRef);

    % Square property
    squareInfo = SquarePropRingSubSubSq(iRef);
    
    % LDFE-center
    quadrature = QDFEcenter(squareInfo);    
    
    % Integration error
    intError = SphInt(quadrature, iRef, intError);
    
    % Angular mesh size
    meshSize(counter) = 1 / sqrt(9 * (iRef + 1) ^ 2 * 24);

    % Print quadrature for PDT
    PrintPDT(quadrature, iRef, fid);
      
    % Update plot counter
    counter = counter + 1;
    
end

% Close PDT quadrature file
fclose(fid);

% generate relative error plots
PlotError(intError, meshSize);