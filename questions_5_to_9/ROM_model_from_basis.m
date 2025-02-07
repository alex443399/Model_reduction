function [A_state_matrix, B_input_matrix] = ROM_model_from_basis(basis_in_matrix_form, ...
    physical_data, geo_data)
    
    resolution = size(basis_in_matrix_form, 2);
    
    hX = physical_data.Lx/(resolution-1);
    hY = physical_data.Ly/(resolution-1);
    % inner_products_between_gradients vs inner_products_with_laplacians
    A_state_matrix_base = inner_products_between_gradients(basis_in_matrix_form, hX, hY);
    B_input_matrix_base = inner_products_with_input_indicators(basis_in_matrix_form, ...
    geo_data, hX, hY, physical_data);

    A_state_matrix = A_state_matrix_base * physical_data.kappa / physical_data.rho / physical_data.c;
    B_input_matrix = B_input_matrix_base / physical_data.rho / physical_data.c;
    
end


%% Functions
function inner_product = Inner_Product(x, y, DeltaX, DeltaY)
    inner_product = sum(x.* y,'all') * DeltaX * DeltaY;
end

function matrix = inner_products_with_laplacians(basis_in_matrix_form, DeltaX, DeltaY)
    R = size(basis_in_matrix_form,1);
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

function input_matrix_pre = inner_products_with_input_indicators(basis_in_matrix_form, ...
    geo_data, hX, hY, physical_data)

    R_reduced_model_order = size(basis_in_matrix_form, 1);
    resolution = size(basis_in_matrix_form, 2);

    input_matrix_pre = zeros(R_reduced_model_order, 2);
    
    Indicator1 = indictator_function_as_mesh(resolution, ...
        geo_data.X1, geo_data.Y1, geo_data.W, physical_data);

    Indicator2 = indictator_function_as_mesh(resolution, ...
        geo_data.X2, geo_data.Y2, geo_data.W, physical_data);
    
    for i_basis_index = 1:R_reduced_model_order
        phi_i = squeeze(basis_in_matrix_form(i_basis_index,:,:));
        for j_u_index = 1:2
            if j_u_index == 1
                input_matrix_pre(i_basis_index,j_u_index) = Inner_Product( ...
                    phi_i, Indicator1, hX, hY);
            end
            if j_u_index == 2
                input_matrix_pre(i_basis_index,j_u_index) = Inner_Product( ...
                    phi_i, Indicator2, hX, hY);
            end
        end
    end
 
end

function matrix = inner_products_between_gradients(basis_in_matrix_form, DeltaX, DeltaY)
    R = size(basis_in_matrix_form, 1);
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


