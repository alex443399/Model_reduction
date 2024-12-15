%% We first write T as a matrix to compute SVD
T_data_matrix = reshape(T, [time_steps, resolution^2])';
disp(['N > M: ', num2str(length(T_data_matrix)>length(T_data_matrix(1,:)))])
%% We compute scd
[Y, Sigma, U] = svd(T_data_matrix);
%% We plot the singular values to observe
singular_values = diag(Sigma);
semilogy(singular_values)
%% Chose amount of singular values based on energy
% This is how much percent of energy we want to keep. I chose 1-1e-6
% because that is roughly the noise level for numerical noise. Based on
% that tolerance we compute R.

Target_energy = 1-1e-6;
percent_of_energy_until = cumsum(singular_values)/sum(singular_values);
R_reduced_model_order = find(percent_of_energy_until>Target_energy,1,"first");
disp(['R: ', num2str(R_reduced_model_order)])
%% Basis functions
% We take the first R basis functions and write them as matrices
basis_in_vector_form = Y(:,1:R_reduced_model_order);
basis_in_matrix_form = permute(reshape( ...
    basis_in_vector_form, [resolution, resolution, R_reduced_model_order]), ...
    [3, 1, 2]);
%% Plot
% We plot the first 6 basis functions
figure;
for i = 1:6
    subplot(2,3,i);
    mesh(x_values, y_values, reshape(basis_in_matrix_form(i,:,:),[resolution,resolution]));
    title(string(i) + 'th basis function')
    xlabel('x')
    ylabel('y')
    zlabel('T')
    hold on
end
hold off