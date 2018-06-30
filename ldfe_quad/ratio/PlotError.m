% Plot integration error with mesh size
function PlotError(intError, meshsize)

% Retrieve information
omegaXerr = intError.omegaXerr;
omegaX2err = intError.omegaX2err;
omegaXomegaYerr = intError.omegaXomegaYerr;
omegaX3err = intError.omegaX3err;
omegaX4err = intError.omegaX4err; 
omegaX3omegaYomegaZerr = intError.omegaX3omegaYomegaZerr ;
omegaX2omegaY2omegaZ2err = intError.omegaX2omegaY2omegaZ2err;

% omega-x integration error
figure
ratio = loglog(meshsize, omegaXerr, 'g-s');
set(ratio, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
grid on
set(gca, 'FontSize', 14)
xlabel('Mesh Size', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
title('Relative Error of \mu', 'FontSize', 18)
legend('LDFE-ratio', 'Location', 'Southeast')

% omega-x^2 integration error
figure
ratio = loglog(meshsize, omegaX2err, 'g-s');
set(ratio, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
grid on
set(gca,'FontSize', 14)
xlabel('Mesh Size', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
title('Relative Error of \mu^2', 'FontSize', 18)
legend('LDFE-ratio', 'Location', 'Southeast')

% omega-x*omega-y integration error
m4 = 4;
b4 = log(omegaXomegaYerr(1));
figure
ratio = loglog(meshsize, omegaXomegaYerr, 'g-s');
set(ratio, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
fourthOrder = loglog(meshsize, exp(m4.*log(meshsize / meshsize(1)) + b4), 'r');
set(fourthOrder, 'LineWidth', 1.25);
grid on
set(gca,'FontSize', 14)
xlabel('Mesh Size', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
title('Relative Error of \mu\eta', 'FontSize', 18)
legend('LDFE-ratio', 'Fourth Order', 'Location', 'Southeast')

% omega-x^3 integration error
m4 = 4;
b4 = log(omegaX3err(1));
figure
ratio = loglog(meshsize, omegaX3err, 'g-s');
set(ratio, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
fourthOrder = loglog(meshsize, exp(m4.*log(meshsize / meshsize(1)) + b4), 'r');
set(fourthOrder, 'LineWidth', 1.25);
grid on
set(gca,'FontSize', 14)
xlabel('Mesh Size', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
title('Relative Error of \mu^3', 'FontSize', 18)
legend('LDFE-ratio', 'Fourth Order', 'Location', 'Southeast')

% omega-x^4 integration error
m4 = 4;
b4 = log(omegaX4err(1));
figure
ratio = loglog(meshsize, omegaX4err, 'g-s');
set(ratio, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
fourthOrder = loglog(meshsize, exp(m4.*log(meshsize / meshsize(1)) + b4), 'r');
set(fourthOrder, 'LineWidth', 1.25);
grid on
set(gca,'FontSize', 14)
xlabel('Mesh Size', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
title('Relative Error of \mu^4', 'FontSize', 18)
legend('LDFE-ratio', 'Fourth Order', 'Location', 'Southeast')

% omega-x^3*omega-y*omega-z integration error
m4 = 4;
b4 = log(omegaX3omegaYomegaZerr(1));
figure
ratio = loglog(meshsize, omegaX3omegaYomegaZerr, 'g-s');
set(ratio, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
fourthOrder = loglog(meshsize, exp(m4.*log(meshsize / meshsize(1)) + b4), 'r');
set(fourthOrder, 'LineWidth', 1.25);
grid on
set(gca,'FontSize', 14)
xlabel('Mesh Size', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
title('Relative Error of \mu^2\eta^2\epsilon', 'FontSize', 18)
legend('LDFE-ratio', 'Fourth Order',...
  'Location', 'Southeast')

% omega-x^2*omega-y^2*omega-z^2 integration error
m4 = 4;
b4 = log(omegaX2omegaY2omegaZ2err(1));
figure
ratio = loglog(meshsize, omegaX2omegaY2omegaZ2err, 'g-s');
set(ratio, 'LineWidth', 1.25, 'MarkerSize', 10, 'color', [0 0.5 0]);
hold on
fourthOrder = loglog(meshsize, exp(m4.*log(meshsize / meshsize(1)) + b4), 'r');
set(fourthOrder, 'LineWidth', 1.25);
hold on
grid on
set(gca,'FontSize', 14)
xlabel('Mesh Size', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
title('Relative Error of \mu^2\eta^2\epsilon^2', 'FontSize', 18)
legend('LDFE-ratio', 'Fourth Order', 'Location', 'Southeast')
end