function plot_mesh_experiment(Lx, Ly, A_result, t_result, simulation_time_in_seconds, resolution_mesh_grids)
    % configuring the data
    size_A_result = size(A_result);
    t_steps = size_A_result(1);
    K_plus_one = size_A_result(2);
    L_plus_one = size_A_result(3);
    % We compute the initial conditions, this is so we can bound the Z axis
    % on the final part. Since if there is no energy into the system, the
    % initial conditions have the max and min because of how the heat
    % equation works
    A0 = reshape(A_result(1,:,:),[K_plus_one,L_plus_one]);
    
    [~, ~, T_first_sample] =...
        convert_back_from_fourier_efficient(A0, Lx, Ly, ...
        resolution_mesh_grids);

    axis_data = [0, Lx, 0, Ly, min(T_first_sample, [], 'all'), max(T_first_sample, [], 'all')];
    
    % Find the best sample
    [~ , closest_to_desired_time] = min(abs(t_result-simulation_time_in_seconds));

    % we chose a sample
    t_sample = closest_to_desired_time;
    A_sample = reshape(A_result(t_sample,:,:),[K_plus_one,L_plus_one]);
    [x_result_sample, y_result_sample, T_result_sample] =...
        convert_back_from_fourier_efficient(A_sample, Lx, Ly, ...
        resolution_mesh_grids);
    
    % And we graph it
    mesh(x_result_sample, y_result_sample, T_result_sample')
    axis(axis_data)
    title('Simulation at t=' + string(t_result(t_sample))+ 's')
    xlabel('x');
    ylabel('y');
    zlabel('T');
end