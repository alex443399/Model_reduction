L_x = 0.2;
L_y = 0.3;
l_x = 0.04;
l_y = 0.08;
K = 16;
L = 16;
kappa = 148;


Lengths = [L_x, L_y];
X_center = [L_x/4, L_x/4*3];
Y_center = [L_y/2, L_y/2];
W = 0.05;

e_ij = e_mat(K,L, L_x, L_y, l_x, l_y) * (2300*680 - 2328*700);
e_00 = ((2300*680 - 2328*700) * l_x*l_y + 2328*700 * L_x*L_y)/sqrt(L_x*L_y);


T_nm_kl = rank_4_coefficient_tensor(K,L, e_ij, L_x, L_y);
% T_nm_kl(1,1,1,1) = T_00^00
% T_nm_kl(n+1,m+1,k+1,l+1) = T_kl^nm
%% Convert rank-4 tensor to matrix
% let p~nm, q~kl, then we want T_rank_2 to map from R^KL -> R^KL
% Where (T_rank_2)_pq = (T_rank_4)^nm_kl
vector_size = (K+1)*(L+1);
T_matrix = zeros(vector_size, vector_size);
A_matrix = zeros(vector_size, vector_size);
for p = 1:vector_size
    for q = 1:vector_size
        % p = (K+1) l + k + 1
        % q = (K+1) m + n + 1
        k = mod(p-1,K+1);
        l = (p-1-k)/(K+1);

        n = mod(q-1, K+1);
        m = (q-1-n)/(K+1);

        T_matrix(q,p) = T_nm_kl(n+1,m+1,k+1,l+1); % The zero occupies a space
        if p == q
            A_matrix(q,p) = - kappa * ((n/L_x)^2 + (m/L_y)^2);
        end
    end
end

B_matrix = zeros(vector_size, 2);
for q = 1:vector_size
    n = mod(q-1, K+1);
    m = (q-1-n)/(K+1);
    
    omega_x = n*pi/L_x;
    omega_y = m*pi/L_y;
    
    for i = 1:2
         x_term = sqrt(8*L_x)/n/pi * cos(omega_x * X_center(i)) * sin(omega_x * W/2);
         if n == 0
             x_term = W / sqrt(L_x);
         end
         y_term = sqrt(8*L_y)/m/pi * cos(omega_y * Y_center(i)) * sin(omega_y * W/2);
         if m == 0
             y_term = W / sqrt(L_y);
         end
         B_matrix(q,i) = x_term * y_term;
    end
end
E_matrix = eye(vector_size) * e_00 + T_matrix;
sys_weighted = dss(A_matrix, B_matrix, eye(vector_size), 0, E_matrix);
%%
save('ss_with_weigths', "sys_weighted");
%% Simulating
A0 = cubic_initial_conditions(K,L,L_x,L_y);
a0_p = reshape(A0, [vector_size, 1]);

t_span = 0:0.1:600;
u = 1e5 * [1;-1] * sin(t_span/50);
%%
[a_p_sim, ~] = lsim(sys_weighted, u, t_span, a0_p);
%% Reshaping
A_sim = reshape(a_p_sim, [length(t_span),K+1,L+1]);
%% Evaluating
t_show = 1000;
A_t = squeeze(A_sim(t_show,:,:));

[x_values, y_values, T_t] = convert_back_from_fourier_efficient(A_t,L_x,L_y,500);

mesh(x_values, y_values, T_t')
xlabel('x')
ylabel('y')
%% Better plotting

Delta_t = 1;

figure;
for i = 1:6
    subplot(2,3,i);
    plot_mesh_experiment(L_x, L_y, A_sim, ...
    t_span, Delta_t*(i-1), 100)
    hold on
end
hold off
%% Checking e
full_e = [
    e_00 zeros(1,L);
    zeros(L,1) e_ij];
[~, ~, e_in_mesh] = convert_back_from_fourier_efficient(full_e,L_x,L_y,500);
mesh(e_in_mesh);
%% Functions

function e_vec = half_e(K, L_space, l)
    i_vec = 1:K;
    angle_vec = i_vec * pi / 2;
    e_vec = cos(angle_vec) .* sin(angle_vec * l/L_space) ./ angle_vec * sqrt(2* L_space);
end

function e_mat = e_mat(K,L, L_x, L_y, l_x, l_y)
    e_mat = half_e(K,L_x,l_x)' * half_e(L,L_y,l_y);
end

function t = rank_4_coefficient_tensor(K,L, e_ij, L_x, L_y)
    t = zeros(K+1,L+1,K+1,L+1); % First two for n,m last two for k,l
    for n=1:K+1
        for m=1:L+1
            for k=1:K+1
                for l=1:L+1
                    t(n,m,k,l) = Term(k-1,l-1,n-1,m-1, L_x, L_y, K, L, e_ij);
                end
            end
        end
    end
end

function t = Term(k,l,n,m, L_x, L_y, K, L, e_ij)
    t = 0;
    for i = 1:K
        for j = 1:L
            t = t + triple_product_2D(i,j,k,l,n,m,L_x,L_y) * e_ij(i,j);
        end
    end
end

function p = triple_product_2D(i,j,k,l,n,m,L_x,L_y)
    p = triple_product_1D(i,k,n, L_x) * triple_product_1D(j,l,m, L_y);
end

function p = triple_product_1D(i,k,n, L)
    assert (i >= 1)
    assert (k>=0)
    assert (n>=0)
    factor = 1/sqrt(L);
    if k > 0 && n >0
        if n == k+i
            p = factor/sqrt(2);
        elseif n== abs(k-i)
            p = factor/sqrt(2);
        else
            p = 0;
        end
    elseif k>0 && n==0
        if k == i
            p = factor;
        else
            p=0;
        end
    else %k == 0
        if n == i
            p = factor;
        else
            p=0;
        end
    end
end
