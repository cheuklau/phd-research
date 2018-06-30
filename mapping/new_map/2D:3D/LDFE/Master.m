%% Preliminary
restoredefaultpath;
%clear;
clc;
close force all;

%% Use fixup
% 0 = no
% 1 = yes
useFixUp = 1;

%% Linear program selection
% 1 = Revised Simplex
% 2 = Bartels-Golub
% 3 = Forrest-Tomlin
% 4 = Interior Point
LPtype = 1;

%% LDFE quadrature orders
orderFrom = 4; % mapping from
orderTo   = 3; % mapping to

%% Generate LDFE quadratures
quadFrom = GenQuad(orderFrom);
quadTo   = GenQuad(orderTo);

%% Test function
% 1 = Linear function
% 2 = 5th order function
% 3 = Spherical source function
% 4 = Delta-gamma delta-theta patch
% 5 = Delta-gamma delta-theta patch with step
functionNumber = 4;

%% Test -- alpha factor
factor_per = 0.0;

%% Spherical source function information (if necessary)
problem.source_radius   = 1.0; % Source radius
problem.source_strength = 1.0; % Source strength
problem.source_sigma_t  = 1.0; % Source total cross section
problem.detector_pos    = [2, 2, 2]; % Particle position

%% Delta-gamma delta-theta information (if necessary)
problem.min_gamma = pi / 8;
problem.max_gamma = 3 * pi / 8;
problem.min_theta = pi/ 8;
problem.max_theta = 3 * pi / 8;

%% Analytic solution for quadrature you are mapping from
[solFrom, limits] = GenSol(quadFrom, functionNumber, problem, factor_per);

%% Mapped solution for quadrature you are mapping to
solTo = Map(solFrom, limits, orderFrom, orderTo, useFixUp, LPtype);

%% Plot mapped solution
PlotSol(quadTo, solTo, functionNumber, orderFrom, orderTo, problem, limits, useFixUp, factor_per);

%% Verify 0th and 1st moments are preserved
CheckMoments(quadFrom, solFrom, quadTo, solTo);