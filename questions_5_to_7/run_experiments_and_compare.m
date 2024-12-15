function [comparisons, relative_errors] = run_experiments_and_compare( ...
    physical_data,geo_data,experiment_data,K_list,L_list)
    % given Ks and Ls, the comparisons are made with respect to the last
    % ones, which we assume are the largest K, L
    
    
    t_max = 60*10;% 10 minutes
    number_exp = length(K_list); % this is the number of experiments
    % We verify that everything is in order
    verify_lists(K_list, L_list)
    
    % First round of experiments to determine dt
    DT = find_DT(physical_data,geo_data,experiment_data,K_list,L_list, t_max);
    number_of_steps = ceil(t_max/DT)+1;
    disp('time discretization completed')

    % Ex: if Dt=0.6 and we need to go from 0 to 1, we need two samples.
    % Plus, we need the sample at t=0.
    t_span_new = linspace(0, t_max, number_of_steps);

    % We now compute the results for the last experiment, so we can compare
    disp('Computing results for (K,L) = (' + string(K_list(end)) + ',' + string(L_list(end)) + ')');
    [A_result_last, ~] = run_experiment( ...
        physical_data,geo_data,experiment_data,K_list(end),L_list(end),t_span_new);

    % Absolute error as a function of model and time step
    comparisons = zeros(number_exp-1, number_of_steps);
    % Norm of the system so we can normalize later to obtain relative error
    norms_for_normalizing = sqrt(sum(A_result_last.^2,'all'));
    
    for expe = 1:number_exp-1
        K_expe = K_list(expe);
        L_expe = L_list(expe);
        
        % Re run the experiment with uniform time step
        [A_result_expe, ~] = run_experiment( ...
        physical_data,geo_data,experiment_data,K_expe,L_expe,t_span_new);
        
        % Find error for each time step
        for t_step = 1:number_of_steps
            % Reshape the tensors to matrices
            large_A = reshape(A_result_last(t_step,:,:),[K_list(end)+1,L_list(end)+1]);
            current_A = reshape(A_result_expe(t_step,:,:),[K_expe+1,L_expe+1]);
            
            % Get the error and the norm
            error_expe_t = sqrt(get_error_sq_from_matrices(large_A, current_A, K_expe, L_expe));
            comparisons(expe, t_step) = error_expe_t;
        end
        disp('Comparing with (K,L) = (' + string(K_list(expe)) + ',' + string(L_list(expe)) + ')')
    end
    % Normalizing
    relative_errors = comparisons/norms_for_normalizing;

end

function verify_lists(K_list, L_list)
    assert(length(K_list) == length(L_list), 'Amount of K and L should be the same');
    assert(K_list(end) == max(K_list), 'The last K should be the largest');
    assert(L_list(end) == max(L_list), 'The last L should be the largest');
end

function DT = find_DT(physical_data,geo_data,experiment_data,K_list,L_list, t_max)
    t_span = [0 t_max];
    number_exp = length(K_list);

    dt_min_for_exp = zeros(number_exp,1);
    for expe = 1:number_exp
        % We determine minimum timestep for each experiment
        K = K_list(expe);
        L = L_list(expe);
        disp('Computing minimum dt for (K,L) = (' + string(K) + ',' + string(L) + ')')
        % We run the experiment
        [~, t_result_expe] = run_experiment(physical_data,geo_data, experiment_data, K, L, t_span);
        
        % Find the minimum time step for the experiment
        dt_min_for_exp(expe) = min(diff(t_result_expe));
    end
    % We determine minimum timestep across all experiments
    DT = min(dt_min_for_exp);
end

function error_sq = get_error_sq_from_matrices(large_A, current_A, K_expe, L_expe)
    % Since the matrices are of different sizes we cannot substract them
    % We take a submatrix of large_A that we can substract
    submatrix_of_large_matrix = large_A(1:K_expe+1,1:L_expe+1);

    % Computing norms
    % The first norm corresponds to the terms that are in the large model
    % but not in the smaller one. To do this we compute the frobenius norm
    % of the larg matrix and substract the submatrix that overlaps with the
    % other
    higher_order_terms_error = norm(large_A, 'fro')...
        - norm(submatrix_of_large_matrix, 'fro');
    % For the shared terms, we substract them and take the norm
    lower_order_terms_error = norm(submatrix_of_large_matrix-current_A, 'fro');
    % We add the errors
    error_sq = lower_order_terms_error^2 + higher_order_terms_error^2;
    
end