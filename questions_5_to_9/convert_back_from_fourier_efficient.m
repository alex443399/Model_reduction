function [x_values, y_values, T] = convert_back_from_fourier_efficient( ...
    A,Lx,Ly,resolution)
    sizeA_minus_one = size(A)-1;
    K = sizeA_minus_one(1);
    L = sizeA_minus_one(2);

    % points where we will compute the mesh
    x_values = ((1:resolution)-1)/(resolution-1)*Lx;
    y_values = ((1:resolution)-1)/(resolution-1)*Ly;
    
    Dx = x_values(2)-x_values(1);
    Dy = y_values(2)-y_values(1);
    
    % Where we store the mesh
    T = zeros(resolution, resolution);

    for k = 0:K
        for l = 0:L
            % This is \phi_{kl} evaluated in every sample point in
            % x_values, y_values
            summand = eval_basis(k, x_values, Lx)' * eval_basis(l, y_values, Ly);
            % We weight it with its fourier coefficient and add it
            T = T + A(k+1,l+1) * summand;
        end
    end
end