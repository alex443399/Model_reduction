% % We store the variables together for convenience
% physical_data.Lx = 0.2;
% physical_data.Ly = 0.3;
% physical_data.rho = 2328; % 2300;
% physical_data.c = 700; % 680
% physical_data.kappa = 148;
% physical_data.Tamb = 309;
% 
% geo_data.X1 = physical_data.Lx/4;
% geo_data.Y1 = physical_data.Ly/2;
% geo_data.X2 = 3*physical_data.Lx/4;
% geo_data.Y2 = physical_data.Ly/2;
% geo_data.W = 0.05;
% 
% experiment_data.Initial_conditions = 'Cubic'; % It has to be in ['Jump', 'Cubic', 'Zero]
% experiment_data.u_functions = 'Zero'; % It has to be in ['Sines', 'Zero']
% 
% %%
% t_span_equi = linspace(0,600,1817);
% 
% [a_1, t_1, v_1] = run_experiment(physical_data,geo_data, experiment_data, 2, 2, t_span_equi);
% [a_large, t_large, v_large] = run_experiment(physical_data,geo_data, experiment_data, 16, 16, t_span_equi);
% %%
% error_max_question = a_1(:,6) - a_large(:,20);
% plot(error_max_question)
% 
% %%
% error_tensor = A_1-A_large(:,1:3,1:3);
% plot(max(abs(error_tensor), [], [2,3]));
% %%
% heatmap(reshape(error_tensor(884,:,:),[3,3]))
% plot(error_tensor(:,3,2))
%% MORE IN DEPTH

% a2 = reshape(cubic_initial_conditions(2,2,physical_data.Lx,physical_data.Ly), [9 1]);
% a16 = reshape(cubic_initial_conditions(16,16,physical_data.Lx,physical_data.Ly), [17*17 1]);
% %%
% a2(6) == a16(20)
% v_1(6) == v_large(20)
% %%
% f1 = @(t,a) v_1.* a;
% [~, b_1] = ode45(@(t,a) a.*v_1, linspace(0,600,1817), a2);
% %%
% [~, b_large] = ode45(@(t,a) v_large.*a, linspace(0,600,1817), a16);
% %%
% [~, b_singleton] = ode45(@(t,a) v_1(6) * a, linspace(0,600,1817), a2(6));
% %%
% plot(max(a_1-b_1,[],2))
% %%
% plot(max(a_large-b_large,[],2))
% %%
% plot(b_1(:,6)-b_large(:,20))
% %%
% plot(b_1(:,6)-b_singleton)
%% Isolating the error
t_span = linspace(0,600,1817);

coefs = [-1.4023; -0.7648];
initial_conditions = [8; -12];


[~, b_2] = ode45(@(t,a) coefs.*a, t_span, initial_conditions);
[~, b_1] = ode45(@(t,a) coefs(1) .* a, t_span, initial_conditions(1));

max(abs((b_2(:,1)-b_1)))
