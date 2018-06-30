%% Preliminary
restoredefaultpath;
clear;
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

%% QDFE-SQ quadrature orders
orderFrom = 3; % mapping from
orderTo   = 2; % mapping to

%% Generate QDFE-SQ quadratures
quadFrom = GenQuad(orderFrom);
quadTo   = GenQuad(orderTo);

%% Test function
% 1 = 1st-order function
% 2 = 5th-order function
% 3 = Spherical source function
% 4 = Discontinuous function
functionNumber = 2;

%% Additional problem information
if functionNumber == 3
    problem.source_radius   = 1.0;       % Source radius
    problem.source_strength = 1.0;       % Source strength
    problem.source_sigma_t  = 1.0;       % Source total cross section
    problem.detector_pos    = [2, 2, 2]; % Detector position    
elseif functionNumber == 4
    problem.min_gamma = pi / 8;     % Minimum gamma
    problem.max_gamma = 3 * pi / 8; % Maximum gamma
    problem.min_theta = pi/ 8;      % Miminum theta
    problem.max_theta = 3 * pi / 8; % Maximum theta
else
    problem = [];
end

%% Analytic solution for quadrature you are mapping from
[solFrom, limits] = GenSol(quadFrom, functionNumber, problem);

%% Mapped solution for quadrature you are mapping to
solTo = Map(solFrom, limits, orderFrom, orderTo, useFixUp, LPtype);

%% Plot mapped solution
PlotSol(quadTo, solTo, functionNumber, orderFrom, orderTo, problem, limits, useFixUp);

%% Verify 0th and 1st moments are preserved
CheckMoments(quadFrom, solFrom, quadTo, solTo);