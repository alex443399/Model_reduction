function T = convert_from_fourier_to_data_tensor(A_High_dim, resolution, physical_data)
    t_len = size(A_High_dim, 1);
    K = size(A_High_dim, 2)-1;
    L = size(A_High_dim, 3)-1;
    
    x_values = ((1:resolution)-1)/(resolution-1)*physical_data.Lx;
    y_values = ((1:resolution)-1)/(resolution-1)*physical_data.Ly;

    % We initialize the array to store the values, where the axes are t, x, y
    T = zeros(t_len,resolution,resolution);

    for n = 0:K
        for m = 0:L
            % We iterate over the basis functions
            disp(string(n) + ',' + string(m))
            % We sample the basis function in our mesh, we do this outside the
            % time loop to avoid repeated equations
            basis_evaluated = eval_basis(n, x_values, physical_data.Lx)' * eval_basis(m, y_values, physical_data.Ly);
            for t = 1:t_len
                % We compute the contribution for one time step based on
                % fourier coefficient and add it.
                Contributions = A_High_dim(t,n+1,m+1) * basis_evaluated;
                T(t,:,:) = T(t,:,:) + reshape(Contributions, [1, resolution, resolution]);
            end
        end
    end
end

