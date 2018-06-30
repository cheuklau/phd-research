%% Preliminary
restoredefaultpath;
close force all;
clear;
clc;

%% Problem information
prob.x = 35.0;
prob.y = 95.0;
prob.z = 35.0;
prob.sigma_a_source = 1.0;
prob.sigma_a_box    = 1.0;
prob.sigma_a_duct   = 10^-4;
prob.source         = 1.0;

%% Refinement levels
ref_LDFE   = [0 2 3 8 20 40];
ref_QDFE   = [0 1 2 5 10 20 25];
ref_LS     = [4 8 12 16 20 24];
ref_GLC    = [4 16 32 48 64 96 128];
ref_QR     = [2 4 6 8 10 12];
ref_LDFEST = [0 1 2 3 4 5 6 7];

%% Reference solution
%refSol = getRefSol(prob);
refSol = 3.37981E-5;

%% Calculate LDFE error
[intErrorLDFE, intSolLDFE, meshSizeLDFE] = getLDFEerror(ref_LDFE, refSol, prob);

%% Generate QDFE error
[intErrorQDFE, intSolQDFE, meshSizeQDFE] = getQDFEerror(ref_QDFE, refSol, prob);

%% Generate LS error
[intErrorLS, intSolLS, meshSizeLS] = getLSerror(ref_LS, refSol, prob);

%% Generate GLC error
[intErrorGLC, inSolGLC, meshSizeGLC] = getGLCerror(ref_GLC, refSol, prob);

%% Generate QR error
[intErrorQR, intSolQR, meshSizeQR] = getQRerror(ref_QR, refSol, prob);

%% Generate LDFE-ST error
[intErrorLDFEST, intSolLDFEST, meshSizeLDFEST] = getLDFESTerror(ref_LDFEST, refSol, prob);

%% Create error plots
PlotError(intErrorLDFE, meshSizeLDFE, ...
          intErrorQDFE, meshSizeQDFE, ...
          intErrorLS,   meshSizeLS, ...
          intErrorGLC,  meshSizeGLC, ...
          intErrorQR,   meshSizeQR, ...
          intErrorLDFEST, meshSizeLDFEST);