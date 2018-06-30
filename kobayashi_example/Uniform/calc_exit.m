function pathlength = calc_exit(Vmin, P0, P1, x_bounds, y_bounds, z_bounds, p_enter, i)

% Force pathlength to be zero for corner case
force = 0;

% Initialize normals
nx_neg = [-1  0  0];
ny_neg = [ 0 -1  0];
nz_neg = [ 0  0 -1];

% Check face -x
[I, check] = plane_line_intersect(nx_neg,Vmin,P0,P1);
if check ~= 0 && ...
        I(2) >= y_bounds(1) && ...
        I(2) <= y_bounds(2) && ...
        I(3) >= z_bounds(1) && ...
        I(3) <= z_bounds(2)
    
    % Store exiting
    p_exit = I;
    
else
    
    % Check face -y
    [I, check] = plane_line_intersect(ny_neg,Vmin,P0,P1);
    if check ~= 0 && ...
            I(1) >= x_bounds(1) && ...
            I(1) <= x_bounds(2) && ...
            I(3) >= z_bounds(1) && ...
            I(3) <= z_bounds(2)
        
        % Store exiting
        p_exit = I;
        
    else
        
        % Check face -z
        [I, check] = plane_line_intersect(nz_neg,Vmin,P0,P1);
        if check ~= 0 && ...
                I(1) >= x_bounds(1) && ...
                I(1) <= x_bounds(2) && ...
                I(2) >= y_bounds(1) && ...
                I(2) <= y_bounds(2)
            
            % Store exiting
            p_exit = I;
            
            
        else
            
            sprintf('Hit a potential corner case for direction %i \n', i);
            force = 1;
            
        end
        
    end
    
end

% Calculate path length
if force == 0
    
    pathlength = sqrt((p_enter(1) - p_exit(1))^2 + (p_enter(2) - p_exit(2))^2 + (p_enter(3) - p_exit(3))^2);

else
    
    pathlength = 0;
    
end

end