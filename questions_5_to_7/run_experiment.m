function [A_result, t_result] = run_experiment( ...
    physical_data,geo_data, experiment_data, K, L, t_span)
    
    if strcmp(experiment_data.Initial_conditions,'Cubic')
        A0 = cubic_initial_conditions(K,L,physical_data.Lx,physical_data.Ly);
    elseif strcmp(experiment_data.Initial_conditions, 'Jump')

        A0 = jump_initial_conditions(K,L,physical_data.Lx,physical_data.Ly);
    elseif strcmp(experiment_data.Initial_conditions, 'Zero')

        A0 = zeros(k,L);
    else
        assert(False, 'Unknown initial conditions')
    end

    if strcmp(experiment_data.u_functions,'Sines')
        u1_function = @(t) 0.01 * sin(t/50);
        u2_function = @(t) -0.01 * sin(t/50);
    elseif strcmp(experiment_data.u_functions, 'Zero')
        u1_function = @(t) 0;
        u2_function = @(t) 0;
    else
        assert(False, 'Unknown u function')
    end

    [A_result, t_result] = run_experiment_given_initial_conditions( ...
    physical_data,geo_data,K,L,A0,u1_function,u2_function, t_span);
end

function [A_result, t_result] = run_experiment_given_initial_conditions( ...
    physical_data,geo_data,K,L,A0,u1,u2,t_span)
    % We unpack the data
    Lx = physical_data.Lx;
    Ly = physical_data.Ly;
    rho = physical_data.rho;
    c = physical_data.c;
    kappa = physical_data.kappa;
    Tamb = physical_data.Tamb;


    X1 = geo_data.X1;
    X2 = geo_data.X2;
    Y1 = geo_data.Y1;
    Y2 = geo_data.Y2;
    W = geo_data.W;
    
    
    % Contributions of each input to the fourier coefficients
    u_scaling_constant = physical_data.rho*physical_data.c;
    u_1_contributions_to_A_matrix = contribution_single_input_2D(K,L,X1,Y1,W,Lx,Ly)/u_scaling_constant;
    u_2_contributions_to_A_matrix = contribution_single_input_2D(K,L,X2,Y2,W,Lx,Ly)/u_scaling_constant;
    
    u_1_contributions_to_a_vector = reshape( ...
        u_1_contributions_to_A_matrix, [(K+1)*(L+1),1]);
    
    u_2_contributions_to_a_vector = reshape( ...
        u_2_contributions_to_A_matrix, [(K+1)*(L+1),1]);

    % Coefficients for the diff equation
    Coefficients = ((0:K).^2/Lx^2)' + ((0:L).^2/Ly^2);
    vec_coefficients = reshape(Coefficients, [(K+1)*(L+1),1]);
    vector_of_coefficients = (-kappa * pi^2/rho/c) * vec_coefficients;
    
    % We turn it into a vector to make it compatible with a0, and ode45, we
    % could use a diagonal matrix instead of an element wise multiplication

    % Convert initial conditions to vector for ode45 solver
    a0 = reshape(A0,[(K+1)*(L+1),1]);
    % Run ODE

    [t_result, a_result] = ode45( ...
        @(t,a) heat_evolution(t,a, ...
            vector_of_coefficients, u1, u2, ...
            u_1_contributions_to_a_vector, u_2_contributions_to_a_vector), ...
        t_span, a0);
    A_result = reshape(a_result, [length(a_result(:,1)), K+1,L+1]);
    % A_result is rank 3 tensor
    % A_result(t,n,m) = a_{nm}(t)

end

function dadt = heat_evolution(t,a,vector_of_coefficients, u1,u2,u1_cont,u2_cont)
    autonomous = a.*vector_of_coefficients;
    forced_1 = u1(t) * u1_cont;
    forced_2 = u2(t) * u2_cont;
    dadt = forced_1 + forced_2 + autonomous; %
end
