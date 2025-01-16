%% Calling reduced order model
% setting up the matrix of inner products with laplacians of phi
hX = physical_data.Lx/(resolution-1);
hY = physical_data.Ly/(resolution-1);

B_matrix_of_inner_products_with_laplacian = inner_products_between_gradients(basis_in_matrix_form ...
    , R_reduced_model_order, hX, hY);
% Setting up the matrix for the inputs

C_matrix_of_inner_products_with_input = inner_products_with_input_indicators(basis_in_matrix_form, ...
    resolution, geo_data, x_values, y_values, R_reduced_model_order, hX, hY);

%% initial condition

a0 = get_jump_initial_conditions_in_reduced_basis(basis_in_matrix_form, ...
    R_reduced_model_order, resolution, physical_data, hX, hY, x_values, y_values);
%% Simulation
u1_function = @(t) 0.01 * sin(t/50);
u2_function = @(t) -0.01 * sin(t/50);


f_R = @(t,a) (physical_data.kappa * B_matrix_of_inner_products_with_laplacian * a + ...
    C_matrix_of_inner_products_with_input * [u1_function(t); u2_function(t)])/physical_data.rho / physical_data.c; 

% 0.01 * sin(t/50)
[time_steps_reduced_order_model, A_result_reduced_order_model] = ode45( ...
    @(t,a) f_R(t,a), [0 600], a0);

plot(A_result_reduced_order_model)
%% Plotting
axis_data = [0, physical_data.Lx, 0, physical_data.Ly, 0, 1];
desired_times_in_seconds= [0, 10, 20, 30, 40, 50];
figure;
for i = 1:length(desired_times_in_seconds)
    subplot(2,3,i);
    [~ , t_index] = min(abs(time_steps_reduced_order_model-desired_times_in_seconds(i)));
    T_at_time = reduced_basis_to_spatial( ...
        A_result_reduced_order_model(t_index, :), basis_in_matrix_form, resolution);
    
    
    mesh(x_values, y_values, T_at_time');
    title('The temperature distribution at time ' + string(time_steps_reduced_order_model(t_index)) + 's')
    xlabel('x')
    ylabel('y')
    zlabel('T')
    axis(axis_data)
    hold on
end
hold off
%%
plot(A_result_reduced_order_model)
%% implement reduced order model
function dadt = a_coefficients(t, a, k, c, p, R, u)
    dadt = zeros(R,1);
    for i = 1:R
        for j = 1:R
            sum_of_inner_products = sum_of_inner_products + Inner_product(lapl_phi_vector(:,j), basis_in_vector_form(:,i))*a(i);
        end
        dadt(i) = k/(c*p)* sum_of_inner_products + 1/(c*p)*Inner_Product(u,basis_in_vector_form(:,i));
    end
end

function inner_product = Inner_Product(x, y, DeltaX, DeltaY)
    inner_product = sum(x.* y,'all') * DeltaX * DeltaY;
end

function matrix = inner_products_with_laplacians(basis_in_matrix_form, R, DeltaX, DeltaY)
    matrix = zeros(R);
    for i = 1:R
        phi_i = squeeze(basis_in_matrix_form(i,:,:));
        for j = 1:R
            phi_j = squeeze(basis_in_matrix_form(j,:,:));
            % We compute the laplacian
            lapl_phi_j = del2(phi_j, DeltaX, DeltaY);
            % Inner product
            matrix(i,j) = Inner_Product(phi_i, lapl_phi_j, DeltaX, DeltaY);
        end
    end
end

function matrix = inner_products_between_gradients(basis_in_matrix_form, R, DeltaX, DeltaY)
    matrix = zeros(R);
    for i = 1:R
        phi_i = squeeze(basis_in_matrix_form(i,:,:));
        [grad_ix, grad_iy] = gradient(phi_i, DeltaX, DeltaY);
        for j = 1:R
            phi_j = squeeze(basis_in_matrix_form(j,:,:));
            [grad_jx, grad_jy] = gradient(phi_j, DeltaX, DeltaY);
            % We compute the inner products
            matrix(i,j) = -sum(grad_ix.*grad_jx + grad_iy.* grad_jy, 'all') * DeltaX * DeltaY;
        end
    end
end

function matrix = inner_products_with_input_indicators(basis_in_matrix_form, ...
    resolution, geo_data, x_values, y_values, R_reduced_model_order, hX, hY)
    matrix = zeros(R_reduced_model_order, 2);
    
    Indicator1 = zeros(resolution);
    Indicator2 = zeros(resolution);
    
    Center1 = [geo_data.X1; geo_data.Y1];
    Center2 = [geo_data.X2; geo_data.Y2];
    
    for i = 1:resolution
        for j = 1:resolution
            current_postion = [x_values(i); y_values(j)];
            dist1 = norm(Center1 - current_postion, Inf);
            dist2 = norm(Center2 - current_postion, Inf);
            if dist1 <= geo_data.W/2
                Indicator1(i,j) = 1;
            end
            if dist2 <= geo_data.W/2
                Indicator2(i,j) = 1;
            end
        end
    end
    
    for i_basis_index = 1:R_reduced_model_order
        for j_u_index = 1:R_reduced_model_order
            if j_u_index == 1
                matrix(i_basis_index,j_u_index) = Inner_Product( ...
                    squeeze(basis_in_matrix_form(i_basis_index,:,:)), Indicator1, hX, hY);
            end
            if j_u_index == 2
                matrix(i_basis_index,j_u_index) = Inner_Product( ...
                    squeeze(basis_in_matrix_form(i_basis_index,:,:)), Indicator2, hX, hY);
            end
        end
    end
 
end

function a0 = get_jump_initial_conditions_in_reduced_basis(basis_in_matrix_form, ...
    R_reduced_model_order, resolution, physical_data, hX, hY, x_values, y_values)
    a0 = zeros(R_reduced_model_order, 1);
    
    T_R_initial = zeros(resolution);
    center = [physical_data.Lx/2; physical_data.Ly/2];
    W_tilde = 0.5 * min(physical_data.Lx, physical_data.Ly);
    
    for i = 1:resolution
        for j = 1:resolution
            current_postion = [x_values(i); y_values(j)];
            distance = norm(current_postion - center, Inf);
            if distance <= W_tilde/2
                T_R_initial(i,j) = 1;
            end
        end
    end
    
    for i = 1:R_reduced_model_order
        a0(i) = Inner_Product(T_R_initial, ...
            squeeze(basis_in_matrix_form(i,:,:)), hX, hY);
    end
end

function T_profile = reduced_basis_to_spatial(a, basis_in_matrix_form, resolution)
    T_profile = zeros(resolution);
    R = length(a);
    for r = 1:R
        T_profile = T_profile + a(r) * squeeze(basis_in_matrix_form(r,:,:));
    end
end