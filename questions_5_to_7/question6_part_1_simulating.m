% We store the variables together for convenience
physical_data.Lx = 0.2;
physical_data.Ly = 0.3;
physical_data.rho = 2328; % 2300;
physical_data.c = 700; % 680
physical_data.kappa = 148;
physical_data.Tamb = 309;

geo_data.X1 = physical_data.Lx/4;
geo_data.Y1 = physical_data.Ly/2;
geo_data.X2 = 3*physical_data.Lx/4;
geo_data.Y2 = physical_data.Ly/2;
geo_data.W = 0.05;

Initial_conditions = 'Jump';
Input_functions = 'Sines';

experiment_data.Initial_conditions = Initial_conditions; % It has to be in ['Jump', 'Cubic', 'Zero]
experiment_data.u_functions = Input_functions; % It has to be in ['Sines', 'Zero']
%% Run experiment
% Chosing model size 
K = 10;
L = 10;

t_span = [0 60*10]; % 10 mins = 10*60s
[A_result, t_result] = run_experiment(physical_data,geo_data, experiment_data, K, L, t_span);
%% Plot results
% We plot the resulting function at many points in time
graph_for_time = @(t) plot_mesh_experiment(physical_data.Lx, physical_data.Ly, A_result, ...
    t_result, t, 100);
Delta_t = 10;

figure;
for i = 1:6
    subplot(2,3,i);
    graph_for_time(Delta_t*(i-1));
    hold on
end
hold off

%% Compare experiments streamlined
model_sizes = [2,4,8,16,32,48];
aggregating_experiments(physical_data, geo_data, ["Cubic","Jump"], ["Sines", "Zero"], model_sizes);

%% UNUSED
%% Compare experiments based on size
% For the method they both lists must have the same size, and the last
% values must be the largest. We compared the error of all simulation
% results against the largest one.
% We measure the error using the fourier coefficients, rather than the
% spatial domain.
K_list = [8,48]; % [2,4,8,16,32,48]
L_list = [8,48];

% We obtain absolute errors and the total energy of the signal.
% The absolute error is the energy of the error as a function of time. 
% This works because all time_steps are the same.
[error_abs, norm_large_signal] = run_experiments_and_compare(physical_data, geo_data, experiment_data, K_list, L_list);
%% Absolute error as a funtion of time and size
% We plot the relative error between the smaller (K,L)s and the last one.
% This is over time so we can see how time and system dynamics affect the error 
clf
for i = 1:length(K_list)-1
    name = '(K,L) = (' + string(K_list(i)) + ',' + string(L_list(i)) ...
        + ") vs (K',L') = ("+ string(K_list(end)) + ',' + string(L_list(end)) + ')';
    semilogy(error_abs(i,:),'DisplayName',name)
    hold on
end
legend()
title('Absolute error for T0: ' + ...
    string(Initial_conditions) + ' and u: '+ string(Input_functions))
xlabel('Time step')
ylabel('error')
hold off
%% Aggregating errors
error_total_abs = sqrt(sum(error_abs.^2, 2));
error_total_rel = error_total_abs/norm_large_signal;
file_name = 'Relative errors for T0= ' + string(Initial_conditions) + ...
    ' and u= '+ string(Input_functions) + '.mat';
save(file_name, 'error_rel')
%%
semilogy(K_list(1:end-1),error_total_rel)
%%
% We aggregate the errors for measuring and comparing
disp('log10 of the rms relative error')
log_of_rms_relative_errors = log10(rms(error_rel,2));
for i = 1:length(K_list)-1
    disp(string(K_list(i)) + ' vs ' + string(K_list(end)) + ': ' + string(log_of_rms_relative_errors(i)))
end

%% Functions
function aggregating_experiments(physical_data, geo_data, initial_conditions, input_functions, sizes)
    for initial_condition_index = 1:length(initial_conditions)
        for input_index = 1:length(input_functions)
            % Unpacking everything
            initial_condition = initial_conditions(initial_condition_index)
            input_function = input_functions(input_index)
            experiment_data.Initial_conditions = initial_condition;
            experiment_data.u_functions = input_function;
            % Running_experiment
            [error_abs, norm_large_signal] = run_experiments_and_compare(physical_data, ...
                geo_data, experiment_data, sizes, sizes);
            % Normalizing by timestep
            amount_of_time_steps = length(error_abs(1,:));
            tfinal = 60*10;
            time_step_normalizer = sqrt(tfinal/amount_of_time_steps);

            % Aggregating
            error_total_abs = sqrt(sum(error_abs.^2, 2)) * time_step_normalizer;
            error_total_rel = error_total_abs/(norm_large_signal*time_step_normalizer);

            % Storing
            file_name_suffix = ' errors, T0= ' + string(initial_condition) + ...
            ' and u= '+ string(input_function) + '.mat';
            save('abs' + file_name_suffix, 'error_total_abs');
            save('rel' + file_name_suffix, 'error_total_rel');
        end
    end
end