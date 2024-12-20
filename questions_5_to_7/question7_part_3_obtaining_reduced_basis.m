%% We first write T as a matrix to compute SVD
T_data_matrix = reshape(T, [time_steps, resolution^2])';
disp(['N > M: ', num2str(length(T_data_matrix)>length(T_data_matrix(1,:)))])
%% We compute scd
[Y, Sigma, U] = svd(T_data_matrix);
%% We plot the singular values to observe
singular_values = diag(Sigma);
loglog(singular_values)
title("Singular values")
xlabel("i")
ylabel("\sigma_i")
%% Chose amount of singular values based on energy
% This is how much percent of energy we want to keep. I chose 1-1e-6
% because that is roughly the noise level for numerical noise. Based on
% that tolerance we compute R.

Target_energy = 1e-2;
cummulative_energy = cumsum(singular_values.^2);
percent_of_energy_missing = 1 - cummulative_energy/cummulative_energy(end);

semilogy(percent_of_energy_missing)
xlabel("r")
ylabel("M(r)")
title("Fraction of missing energy")

%%
R_reduced_model_order = find(percent_of_energy_missing<Target_energy,1,"first");
disp(['R: ', num2str(R_reduced_model_order)])
%% Basis functions
% We take the first R basis functions and write them as matrices
basis_in_vector_form = Y(:,1:R_reduced_model_order);
basis_in_matrix_form = permute(reshape( ...
    basis_in_vector_form, [resolution, resolution, R_reduced_model_order]), ...
    [3, 1, 2]);
%% Plot
% We plot the first 6 basis functions
clf
figure;
for i = 1:6
    subplot(2,3,i);
    mesh(x_values, y_values, reshape(basis_in_matrix_form(i,:,:),[resolution,resolution]));
    title('Basis function ' + string(i))
    xlabel('x')
    ylabel('y')
    zlabel('T')
    hold on
end
hold off