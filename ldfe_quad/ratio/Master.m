% Master file for cube-based LDFE quadrature sets
close force all;
clear;
clc;

% Add path the spherical harmonic functions
path(path, './SH Integration');

% Refinement level
nRef = 10;

% Initialize relative error
intError = [];

% Write to data file
fid = fopen('LDFEQUADRATIO.DAT', 'at');

% Go through each refinement level
counter = 1;

for iRef = 0 : 1 : nRef
    
    fprintf('Currently on refinement level: %i \n', iRef);
    
    if (iRef == 0)
        
        % Square property
        squareInfo = SquareProp(iRef);
        
    else
        
        % Square property
        squareInfo = SquarePropRingSubSubSq(iRef);
        
    end
    
    % LDFE-ratio
    quadrature = LDFEratioBisection(squareInfo, 4);
    
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