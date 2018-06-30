function plot_tri(tri)

%% Plot sub-square lines

% Sub-square information
length  = 1/ sqrt(3);
x_min   = 0;
x_max   = length;
y_min   = 0;
y_max   = length;
counter = 1;

% Top division
x = linspace(x_min, x_max);
y = y_max;
gamma = atan(x / length);
theta = (pi / 2) - atan(y ./ (length * sec(gamma)));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face two
gamma = (pi / 2) - atan(y ./ length);
theta = (pi / 2) - atan(sqrt(3) .* x ./ sqrt(1 + 3 .* y .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face three
gamma = atan(y ./ x);
theta = (pi / 2) - atan(length ./ sqrt(y .^ 2 + x .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;

% Right division
x = x_max;
y = linspace(y_min, y_max);
gamma = atan(x / length);
theta = (pi / 2) - atan(y ./ (length * sec(gamma)));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face two
gamma = (pi / 2) - atan(y ./ length);
theta = (pi / 2) - atan(sqrt(3) .* x ./ sqrt(1 + 3 .* y .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face three
gamma = atan(y ./ x);
theta = (pi / 2) - atan(length ./ sqrt(y .^ 2 + x .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;

% Bottom division
x = linspace(x_min, x_max);
y = y_min;
gamma = atan(x / length);
theta = (pi / 2) - atan(y ./ (length * sec(gamma)));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face two
gamma = (pi / 2) - atan(y ./ length);
theta = (pi / 2) - atan(sqrt(3) .* x ./ sqrt(1 + 3 .* y .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face three
gamma = atan(y ./ x);
theta = (pi / 2) - atan(length ./ sqrt(y .^ 2 + x .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;

% Left division
x = x_min;
y = linspace(y_min, y_max);
gamma = atan(x / length);
theta = (pi / 2) - atan(y ./ (length * sec(gamma)));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face two
gamma = (pi / 2) - atan(y ./ length);
theta = (pi / 2) - atan(sqrt(3) .* x ./ sqrt(1 + 3 .* y .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);
counter = counter + 1;
% Rotate to face three
gamma = atan(y ./ x);
theta = (pi / 2) - atan(length ./ sqrt(y .^ 2 + x .^ 2));
omegaXsub{counter} = cos(gamma) .* sin(theta);
omegaYsub{counter} = sin(gamma) .* sin(theta);
omegaZsub{counter} = cos(theta);

% Plot each line
for k = 1 : counter    
    plot3(omegaXsub{k}, omegaYsub{k}, omegaZsub{k}, 'k', 'LineWidth', 1.25)    
    hold on    
end

%% Plot triangular integration patch

% Triangular patch information
v_top = tri.v_top; % Top of triangle (fraction of base triangle height)
v_bot = tri.v_bot; % Bottom of triangle (fraction of base triangle height)
u_tip = tri.u_tip; % Position of triangle tip (fraction of base at tip)
flip  = tri.flip;  % Triangle orientation

% Base triangle parameters
h = sqrt(6) / 2;     % Height
b = 2 * h / sqrt(3); % Base length (different than dissertation)

% Adjust triangle information based on height and base length
v_top       = v_top * h; % Top of triangle
v_bot       = v_bot * h; % Bottom of triangle
u_tip_tilde = u_tip * (b / 2); % u-tilde for tip

% Gauss order
N = 64;

% Divide v using N-point Gauss quadrature
[fv, ~] = lgwt(N, v_bot, v_top);

% Go over each v quadrature point
counter = 1;
for i = 1 : N
    
    % Determine u bounds for rightside-up triangle
    if flip == 0        
        u_min = (fv(i) - (v_top - (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));        
        u_max = (-1 * fv(i) + (v_top + (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));
        
        % Determine u bounds for upside-down triangle
    else       
        u_min = (-1 * fv(i) + (v_bot + (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));        
        u_max = (fv(i) - (v_bot - (2 * h / b) * u_tip_tilde)) / (h * (1 - fv(i) / h));        
    end
    
    % Divide u using N-point Gauss quadrature
    [fu, ~] = lgwt(N, u_min, u_max);
    
    % Go through each u quadrature point
    for j = 1 : N
        
        % Calculate solid angles based on local triangular coordinates
        gamma = (pi / 4) - atan(-fu(j));
        theta = acos(sqrt(2) * fv(i) / sqrt(h ^ 2 * (1 + fu(j) ^ 2) - ...
            2 * h * (1 + fu(j) ^ 2) * fv(i) + (3 + fu(j) ^ 2) * fv(i) ^ 2));
        
        % Calculate directional cosines
        omegaX(counter) = cos(gamma) * sin(theta);
        omegaY(counter) = sin(gamma) * sin(theta);
        omegaZ(counter) = cos(theta);        
        counter = counter + 1;
        
    end
    
end

% Plot triangular integration patch
scatter3(omegaX, omegaY, omegaZ, 10, 'r', 'filled');
set(gca,'FontSize', 14)
xlabel('x', 'FontSize', 18)
ylabel('y', 'FontSize', 18)
zlabel('z', 'FontSize', 18)

end