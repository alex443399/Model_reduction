R = 3;
%%
experiment_data.Initial_conditions = 'Jump'; % It has to be in ['Jump', 'Cubic', 'Zero]
experiment_data.u_functions = 'Sines'; % It has to be in ['Sines', 'Zero']

% Load basis
basis_in_matrix_form = load('full_100_v3_basis_in_matrix_form.mat').full_basis_in_matrix_form;
basis_in_matrix_form = basis_in_matrix_form(1:R, :, :);
% Variables
resolution = size(basis_in_matrix_form, 2);

%% Data
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
%% Simulation HDM
K = 16;
L = 16;
[A_High_dim, t_High_dim] = run_experiment(physical_data,geo_data, experiment_data, K, L, [0 60*10]);
%% Convert HDM
Jumps = 50;
HDM_TEMP = get_data_tensor(A_High_dim, t_High_dim, physical_data, ...
    resolution, Jumps);

%% Calling reduced order model
% setting up the matrix of inner products with laplacians of phi
[ROM_state_matrix, ROM_input_matrix] = ROM_model_from_basis(basis_in_matrix_form, ...
    physical_data, geo_data);

a0 = reduced_initial_conditions(basis_in_matrix_form, physical_data, experiment_data);
% HDM_MESH_T0 = squeeze(HDM_TEMP(1,:,:));
% a0_THDM = Initial_conditions_given_mesh(basis_in_matrix_form, HDM_MESH_T0, ...
%      physical_data.Lx/(resolution-1), physical_data.Ly/(resolution-1));
u_rom = get_u_function(experiment_data);
%%
% Simulating ROM
[~, a_ROM] = ode45(@(t,a) ROM_state_matrix * a + ROM_input_matrix * u_rom(t), ...
    t_High_dim, a0);
% Mesh
ROM_TEMP = convert_ROM_to_data_tensor(a_ROM, basis_in_matrix_form, Jumps);
%% Comparing meshes
dA = physical_data.Lx * physical_data.Ly /(resolution-1)^2;
DELTA_TEMP = HDM_TEMP - ROM_TEMP;
Error_energy = sum(DELTA_TEMP.^2, [2,3])*dA;
Original_energy = sum(HDM_TEMP.^2, [2,3])*dA;
error_relative = sqrt(sum(Error_energy)/sum(Original_energy))
plot(Error_energy)
%% Plot
figure;
subplot(1,2,1);
mesh(squeeze(HDM_TEMP(2,:,:))');
xlabel('x')
ylabel('y')
title('HDM')
subplot(1,2,2);
mesh(squeeze(ROM_TEMP(2,:,:))');
xlabel('x')
ylabel('y')
title('ROM')
%%
plot(t_High_dim, sum(DELTA_TEMP.^2,[2,3]))

%% Functions
function T_ROM = convert_ROM_to_data_tensor(a_ROM, basis_in_matrix_form, Jumps)
    sampled_time_steps = length(1:Jumps:size(a_ROM,1));
    
    resolution = size(basis_in_matrix_form, 2);
    T_ROM = zeros(sampled_time_steps, resolution, resolution);

    for t_sample = 1:sampled_time_steps
        t_sim = 1 + (t_sample-1) * Jumps;
        a_t = a_ROM(t_sim,:)';
        product = sum(a_t.*basis_in_matrix_form);
        T_ROM(t_sample,:,:) = product;
    end
end
function a0 = Initial_conditions_given_mesh(basis_in_matrix_form, T0, hX, hY)
    R = size(basis_in_matrix_form,1);
    a0 = zeros(R,1);
    for i = 1:R
        phi_i = squeeze(basis_in_matrix_form(i,:,:));
        a0(i) = Inner_Product(T0, phi_i, hX, hY);
    end
end

function inner_product = Inner_Product(x, y, DeltaX, DeltaY)
    inner_product = sum(x.* y,'all') * DeltaX * DeltaY;
end

function c_mesh = cubic_mesh(resolution, physical_data)
    % Sample points
    x_values = ((1:resolution)-1)/(resolution-1)*physical_data.Lx;
    y_values = ((1:resolution)-1)/(resolution-1)*physical_data.Ly;
    
    % Evaluating the function in the sample points
    f = @(x) x.^2.*(0.5-x/3);
    f_eval_x = f(x_values/physical_data.Lx);
    f_eval_y = f(y_values/physical_data.Ly);

    c_mesh = f_eval_x * f_eval_y';
end

function a0 = reduced_initial_conditions(basis_in_matrix_form, physical_data, experiment_data)
    resolution = size(basis_in_matrix_form, 2);
    
    hX = physical_data.Lx/(resolution-1);
    hY = physical_data.Ly/(resolution-1);
    
    if strcmp(experiment_data.Initial_conditions,'Cubic')
        T0 = cubic_mesh(resolution, physical_data);
        a0 = Initial_conditions_given_mesh(basis_in_matrix_form, T0, hX, hY);
    elseif strcmp(experiment_data.Initial_conditions, 'Jump')
        W_init = 0.5 * min(physical_data.Lx, physical_data.Ly);
        T0= indictator_function_as_mesh(resolution, ...
            physical_data.Lx/2, physical_data.Ly/2, ...
        W_init, physical_data);

        a0 = Initial_conditions_given_mesh(basis_in_matrix_form, T0, hX, hY);
    elseif strcmp(experiment_data.Initial_conditions, 'Zero')

        a0 = zeros(size(basis_in_matrix_form,1),1);
    else
        assert(False, 'Unknown initial conditions')
    end
end

function u_fun = get_u_function(experiment_data)
    if strcmp(experiment_data.u_functions,'Sines')
        u_fun = @(t) 1e5 * sin(t/50) * [1;-1];
    elseif strcmp(experiment_data.u_functions, 'Zero')
        u_fun = @(t) [0;0];
    else
        assert(False, 'Unknown u function')
    end
end