function A = cubic_initial_conditions(K,L,Lx,Ly)
    % The cubic equation f as defined below has derivative zero at 0 and 1.
    % Our initial conditions are f(x/Lx)*f(y/Ly)
    % This means that the inner product gets separated into 
    % <f(x/Lx), phi^x> and <f(y/Ly), phi^y> where <.,.> is the 1D L2 inner
    % product. We exploit this for efficiency, computing them independently
    % and then multiplying them out.

    f = @(x) x.^2.*(0.5-x/3);
    resolution =  1e6; % 1e5 seems to be enough for 4 digits after comma
    
    % Sample points
    x_values = ((1:resolution)-1)/(resolution-1)*Lx;
    y_values = ((1:resolution)-1)/(resolution-1)*Ly;
    
    Dx = x_values(2)-x_values(1);
    Dy = y_values(2)-y_values(1);
    % Evaluating the function in the sample points
    f_eval_x = f(x_values/Lx);
    f_eval_y = f(y_values/Ly);
    
    % Storing inner products
    Inner_products_x = zeros(K+1,1);
    Inner_products_y = zeros(L+1,1);
    
    
    for n = 0:K
        cos_x = eval_basis(n, x_values, Lx);
        Inner_products_x(n+1) = Dx * dot(f_eval_x, cos_x); % Riemann summ for inner product
        
    end
    for m = 0:L
        cos_y = eval_basis(m, y_values, Ly);
        Inner_products_y(m+1) = Dy * dot(f_eval_y, cos_y); % Idem
    end

    A = Inner_products_x * Inner_products_y'; % Join results
end

