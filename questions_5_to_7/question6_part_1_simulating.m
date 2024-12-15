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

experiment_data.Initial_conditions = 'Cubic'; % It has to be in ['Jump', 'Cubic', 'Zero]
experiment_data.u_functions = 'Sines'; % It has to be in ['Sines', 'Zero']

% Chosing model size 
K = 10;
L = 10;

%% Run experiment
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

%% Compare experiments based on size
% For the method they both lists must have the same size, and the last
% values must be the largest. We compared the error of all simulation
% results against the largest one.
% We measure the error using the fourier coefficients, rather than the
% spatial domain.
K_list = [2,4,8,16,32, 48]; % [2,4,8,16,32]
L_list = [2,4,8,16,32, 48];

% We obtain absolute errors and relative errors.
% The absolute error is the energy of the error as a function of time. When
% we compute the relative error we divide by the total energy of the
% function over time. We use the largest function. The calc of the total energy is of the form
% sum(A_result, [], 'all'). This works because all time_steps are the same.
[error_abs, error_rel] = run_experiments_and_compare(physical_data, geo_data, experiment_data, K_list, L_list);
%% Relative error as a funtion of time and size
% We plot the relative error between the smaller (K,L)s and the last one.
% This is over time so we can see how time and system dynamics affect the error 
clf
for i = 1:length(K_list)-1
    name = '(K,L) = (' + string(K_list(i)) + ',' + string(L_list(i)) ...
        + ') vs (K,L) = ('+ string(K_list(end)) + ',' + string(L_list(end)) + ')';
    plot(error_rel(i,:),'DisplayName',name)
    hold on
end
hold off
legend()
title('error based on basis size across time')
hold off
%% Aggregating errors
% We aggregate the errors for measuring and comparing
disp('log10 of the rms relative error')
log_of_rms_relative_errors = log10(rms(error_rel,2));
for i = 1:length(K_list)-1
    disp(string(K_list(i)) + ' vs ' + string(K_list(end)) + ': ' + string(log_of_rms_relative_errors(i)))
end
