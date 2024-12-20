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

experiment_data.Initial_conditions = 'Jump'; % It has to be in ['Jump', 'Cubic', 'Zero]
experiment_data.u_functions = 'Sines'; % It has to be in ['Sines', 'Zero']

% Chosing model size 
K = 16;
L = 16;

%% We run the experiment in high dimensions.
[A_High_dim, t_High_dim] = run_experiment(physical_data,geo_data, experiment_data, K, L, [0 60*10]);

%% UNUSED, too much memory space, ~ 200kx64x64
% %% Run experiment to obtain the step size
% t_span_for_finding_Dt = [0 60*10]; % 10 mins = 10*60s
% [~, t_result_for_dt] = run_experiment(physical_data,geo_data, experiment_data, K, L, t_span_for_finding_Dt);
% %%
% Dt = min(diff(t_result_for_dt));
% number_of_equally_spaced_time_samples = ceil(60*10/Dt)+1;
% t_span_uniform = linspace(0,10*60, number_of_equally_spaced_time_samples);
% [A_High_dim, t_High_dim] = run_experiment(physical_data,geo_data, experiment_data, K, L, t_span_uniform);
% %%