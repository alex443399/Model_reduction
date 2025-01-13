%% Calling reduced order model
% setting up the laplacian of phi
hX = physical_data.Lx/(resolution-1);
hY = physical_data.Ly/(resolution-1);
lapl_phi_vector = []; 
sum_of_inner_products = 0;
for j = 1:Reduced_model_order
    phi = squeeze(basis_in_matrix_form(j,:,:));
    lapl_phi = del2(phi, hX, hY);
    lapl_phi_vector = [lapl_phi_vector, reshape(lapl_phi,[], 1)];
end

% initial condition
a0 = zeros(Reduced_model_order, 1);
f = @(x) x.^2.*(0.5-x/3);
    
% Evaluating the function in the sample points
f_eval_x = f(x_values/Lx);
f_eval_y = f(y_values/Ly);
%T_R_initial =  ...
% u = ...

for i = 1:Reduced_model_order
    a0(i) = Inner_Product(T_R_initial, basis_in_vector_form(i));
end

[time_steps, a] = ode45(@a_coefficients, time_steps, a0, [], physical_data.kappa, physical_data.c, physical_data.rho, Reduced_model_order, u); 

basis_in_matrix_form = permute(reshape(a, [resolution, resolution, R_reduced_model_order]), [3, 1, 2]);

temperature_time_points = [1, 10, 50, 100, 1000];
figure;
for i = 1:length(temperature_time_points)
    subplot(2,3,i);
    mesh(x_values, y_values, reshape(basis_in_matrix_form(i,:,:),[resolution,resolution]));
    title('The temperature distribution at time ' + string(i))
    xlabel('x')
    ylabel('y')
    zlabel('T')
    hold on
end
hold off
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

function inner_product = Inner_Product(x, y)
    inner_product = sum(x.* y,'all');
end