% Master file for cube-based LDFE quadrature sets
close force all;
clear;
clc;

% Add path the spherical harmonic functions
path(path, './SH Integration');

% Refinement level
nRef = 127;

% Initialize relative error
intError = [];

% Write to data file
fid = fopen('LDFESQ.DAT', 'at');

% Go through each refinement level
counter = 1;

for iRef = 127 : 1 : nRef
        
    fprintf('Currently on refinement level: %i \n', iRef);
    
    % Square property
    if (iRef == 0) 
                
        squareInfo = SquareProp(iRef);
        
    else 
       
        squareInfo = SquarePropRingSubSubSq(iRef);
        
    end    
    
    % LDFE-ratio for first guess ratios
    [~, rhoInitAll] = genFirstRatios(squareInfo, 2);

    % LDFE-sa
    quadrature = LDFEsa(squareInfo, rhoInitAll);
    
    % Integration error
    intError = SphInt(quadrature, iRef, intError);
    
    % Angular mesh size
    meshSize(counter) = 1 / sqrt(4 * (iRef + 1) ^ 2 * 24);
    
    % Print quadrature for PDT
    PrintPDT(quadrature, iRef, fid);
      
    % Update plot counter
    counter = counter + 1;
   
end

% Close PDT quadrature file
fclose(fid);

% generate relative error plots
PlotError(intError, meshSize);