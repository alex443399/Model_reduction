function T = get_data_tensor(A_High_dim, t_High_dim, physical_data, ...
    resolution, Jumps)
    simulation_time_steps = length(t_High_dim);
    sampled_time_steps = length(1:Jumps:simulation_time_steps);
    
    K = size(A_High_dim,2) - 1;
    L = size(A_High_dim,3) - 1;
    % This are the sample points
    x_values = ((1:resolution)-1)/(resolution-1)*physical_data.Lx;
    y_values = ((1:resolution)-1)/(resolution-1)*physical_data.Ly;
    
    % We initialize the array to store the values, where the axes are t, x, y
    T = zeros(sampled_time_steps,resolution,resolution);
    
    %% We compute T
    for n = 0:K
        for m = 0:L
            % We iterate over the basis functions
            disp(string(n) + ',' + string(m))
            % We sample the basis function in our mesh, we do this outside the
            % time loop to avoid repeated equations
            basis_evaluated = eval_basis(n, x_values, physical_data.Lx)' * eval_basis(m, y_values, physical_data.Ly);
            for t_sample = 1:sampled_time_steps
                t_simulation = 1 + Jumps * (t_sample-1);
                % We compute the contribution for one time step based on
                % fourier coefficient and add it.
                Contributions = A_High_dim(t_simulation,n+1,m+1) * basis_evaluated;
                T(t_sample,:,:) = T(t_sample,:,:) + reshape(Contributions, [1, resolution, resolution]);
            end
        end
    end

end
