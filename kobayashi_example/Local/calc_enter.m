function p_enter = calc_enter(Vmax, P0, P1, x_bounds, y_bounds, z_bounds)

% Initialize normals
nx_pos = [ 1  0  0];
ny_pos = [ 0  1  0];
nz_pos = [ 0  0  1];

% Check face +x
[I, check] = plane_line_intersect(nx_pos,Vmax,P0,P1);
if check ~= 0 && ...
        I(2) >= y_bounds(1) && ...
        I(2) <= y_bounds(2) && ...
        I(3) >= z_bounds(1) && ...
        I(3) <= z_bounds(2)
    
    % Store entering
    p_enter = I;
    
else
    
    % Check face +y
    [I, check] = plane_line_intersect(ny_pos,Vmax,P0,P1);
    if check ~= 0 && ...
            I(1) >= x_bounds(1) && ...
            I(1) <= x_bounds(2) && ...
            I(3) >= z_bounds(1) && ...
            I(3) <= z_bounds(2)
        
        % Store entering
        p_enter = I;
        
    else
        
        % Check face +z
        [I, check] = plane_line_intersect(nz_pos,Vmax,P0,P1);
        if check ~= 0 && ...
                I(1) >= x_bounds(1) && ...
                I(1) <= x_bounds(2) && ...
                I(2) >= y_bounds(1) && ...
                I(2) <= y_bounds(2)
            
            % Store entering
            p_enter = I;
            
        else
            
            p_enter = -1;
            
        end
        
    end
    
end

end