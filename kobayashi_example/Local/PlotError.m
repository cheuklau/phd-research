% Plot integration error with mesh size
function PlotError(intErrorLDFE,   meshSizeLDFE, ...
                   intErrorQDFE,   meshSizeQDFE, ...
                   intErrorLS,     meshSizeLS, ...
                   intErrorGLC,    meshSizeGLC, ...
                   intErrorQR,     meshSizeQR, ...
                   intErrorLDFEST, meshSizeLDFEST)

% Retrieve information
errLDFE   = intErrorLDFE;
errQDFE   = intErrorQDFE;
errLS     = intErrorLS;
errGLC    = intErrorGLC;
errQR     = intErrorQR;
errLDFEST = intErrorLDFEST;

% omega-x * omega-y
m2 = 2;
b2 = log(errLDFE(1));
figure

ldfePlot = loglog(meshSizeLDFE, errLDFE, 'k-+');
set(ldfePlot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on

qdfePlot = loglog(meshSizeQDFE, errQDFE, 'k-o');
set(qdfePlot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on

lsPlot = loglog(meshSizeLS, errLS, 'm-+');
set(lsPlot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on

glcPlot = loglog(meshSizeGLC, errGLC, 'm-o');
set(glcPlot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on

qrPlot = loglog(meshSizeQR, errQR, 'm-s');
set(qrPlot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on

ldfestPlot = loglog(meshSizeLDFEST, errLDFEST, 'b-s');
set(ldfestPlot, 'LineWidth', 1.25, 'MarkerSize', 10);
hold on

secondorder = loglog(meshSizeLDFE, exp(m2.*log(meshSizeLDFE / meshSizeLDFE(1)) + b2), 'r');
set(secondorder, 'LineWidth', 1.25);

grid on
set(gca,'FontSize', 18)
xlabel('Mesh Length', 'FontSize', 24)
ylabel('|Relative Error|', 'FontSize', 24)
legend('LDFE-SQ', 'QDFE-SQ', 'LS', 'GLC', 'QR', 'LDFE-ST', 'Second-order', 'Location', 'Southeast')

end