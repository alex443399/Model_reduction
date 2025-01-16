% Preprocessing the data
% Resolution is the mesh size in spatial coordinates, for example, if it is
% 100 that means that we sample 100x100 points, 100 on each axis
resolution = 100;

time_steps = length(t_High_dim);

% This are the sample points
x_values = ((1:resolution)-1)/(resolution-1)*physical_data.Lx;
y_values = ((1:resolution)-1)/(resolution-1)*physical_data.Ly;

% We initialize the array to store the values, where the axes are t, x, y
T = zeros(time_steps,resolution,resolution);

%% We compute T
for n = 0:K
    for m = 0:L
        % We iterate over the basis functions
        disp(string(n) + ',' + string(m))
        % We sample the basis function in our mesh, we do this outside the
        % time loop to avoid repeated equations
        basis_evaluated = eval_basis(n, x_values, physical_data.Lx)' * eval_basis(m, y_values, physical_data.Ly);
        for t = 1:time_steps
            % We compute the contribution for one time step based on
            % fourier coefficient and add it.
            Contributions = A_High_dim(t,n+1,m+1) * basis_evaluated;
            T(t,:,:) = T(t,:,:) + reshape(Contributions, [1, resolution, resolution]);
        end
    end
end

%% We plot to verify
mesh(x_values,y_values,squeeze(T(1,:,:))');
xlabel('x')
ylabel('y')