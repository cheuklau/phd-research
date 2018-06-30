% Note: LDFE quadratures are pre-generated to save run time

% Preliminary
restoredefaultpath;
clear;
clc;
close force all;

% LDFE quadrature orders
orderFrom = 2; % mapping from
orderTo   = 3; % mapping to

% Generate LDFE quadratures
quadFrom = GenQuad(orderFrom);
quadTo   = GenQuad(orderTo);

% Test function
% 1 = Linear function
% 2 = Quadratic function
% 3 = Spherical source function
functionNumber = 1;

% Spherical source function information
problem.source_radius   = 0.50; % Source radius
problem.source_strength = 1.0;  % Source strength            
problem.source_sigma_t  = 1.0;  % Source total cross section
problem.detector_pos    = [1.0, 1.0, 1.0]; % Particle position

% Analytic solution for quadrature you are mapping from
solFrom = GenSol(quadFrom, functionNumber, problem);

% Mapping algorithm
% 1 = preserve 0th and 1st moments
% 2 = preserve 1st moment using LDFE basis functions
mapping = 1;

% Mapped solution for quadrature you are mapping to
solTo = Map(solFrom, orderFrom, orderTo, mapping);

% Plot mapped solution
PlotSol(quadTo, solTo, functionNumber, orderFrom, orderTo, problem);

% Verify 0th and 1st moments are preserved
CheckMoments(quadFrom, solFrom, quadTo, solTo);