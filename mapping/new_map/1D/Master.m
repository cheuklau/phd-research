% Preliminary
clear;
clc;
close force all;

% LDFE quadrature orders
orderFrom = 3; % mapping from
orderTo   = 4; % mapping to

% Generate LDFE quadratures
quadFrom = GenQuad(orderFrom);
quadTo   = GenQuad(orderTo);

% Test function
% 1 = Step function with discontinuity within an LDFE region
% 2 = Step function with discontinuity at LDFE region boundary
% 3 = Smooth linear function
% 4 = Inverse exponential function with discontinuity at mu=0
% 5 = Preliminary exam problem
functionNumber = 1;

% Analytic solution for quadrature you are mapping from
[solFrom, limits] = GenSol(quadFrom, functionNumber);

% Mapped solution for quadrature you are mapping to
solTo = Map(solFrom, orderFrom, orderTo, limits);

% Plot mapped solution
PlotSol(quadFrom, solFrom, quadTo, solTo, functionNumber);

% Verify 0th and 1st moments are preserved
CheckMoments(quadFrom, solFrom, quadTo, solTo);