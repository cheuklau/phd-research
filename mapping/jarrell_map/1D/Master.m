% Preliminary
clear;
clc;
close force all;

% LDFE quadrature orders
orderFrom = 3; % mapping from 
orderTo   = 2; % mapping to

% Generate LDFE quadratures
quadFrom = GenQuad(orderFrom);
quadTo   = GenQuad(orderTo);

% Test function
% 1 = Step function with discontinuity within an LDFE region
% 2 = Step function with discontinuity at LDFE region boundary
% 3 = Smooth linear function
% 4 = Inverse exponential function with discontinuity at mu=0
functionNumber = 4;

% Analytic solution for quadrature you are mapping from
solFrom = GenSol(quadFrom, functionNumber);

% Mapping algorithm
% 1 = preserve 0th and 1st moments
% 2 = preserve 1st moment using LDFE basis functions
% 3 = preserve 1st moment using CDFE basis functions
mapping = 2;

% Mapped solution for quadrature you are mapping to
solTo = Map(solFrom, orderFrom, orderTo, mapping);

% Plot mapped solution
PlotSol(quadFrom, solFrom, orderFrom, quadTo, solTo, ...
    orderTo, functionNumber);

% Verify 0th and 1st moments are preserved
CheckMoments(quadFrom, solFrom, quadTo, solTo);