function indicator = indictator_function_as_mesh(resolution, center_X, center_Y, W, physical_data)
    indicator = zeros(resolution);
    Center = [center_X; center_Y];
    
    x_values = ((1:resolution)-1)/(resolution-1)*physical_data.Lx;
    y_values = ((1:resolution)-1)/(resolution-1)*physical_data.Ly;
    
    for i = 1:resolution
        for j = 1:resolution
            current_postion = [x_values(i); y_values(j)];
            dist = norm(Center - current_postion, Inf);
            if dist <= W/2
                indicator(i,j) = 1;
            end
        end
    end
end