% Preprocessing the data
% Resolution is the mesh size in spatial coordinates, for example, if it is
% 100 that means that we sample 100x100 points, 100 on each axis
resolution = 100;
Jumps = 100;
T = get_data_tensor(A_High_dim, t_High_dim, physical_data, ...
    resolution, Jumps);


%% We plot to verify
mesh(squeeze(T(60,:,:))');
xlabel('x')
ylabel('y')

%% Functions
