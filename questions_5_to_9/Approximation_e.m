rho1 = 2328;
rho2 = 2300;
c1 = 700;
c2 = 680;
Ly = 0.3;
Lx = 0.2;
lx = 0.04;
ly = 0.08;
e00 = (1/sqrt(Lx*Ly))*(rho1*c1*(Lx*Ly-lx*ly)+rho2*c2*lx*ly);
K = 48;
L = 48;

resolution = 100;

%This are the sample points
x_values = ((1:resolution)-1)/(resolution-1)*Lx;
y_values = ((1:resolution)-1)/(resolution-1)*Ly;

e = zeros(K,L);
phi = zeros(K,L, resolution, resolution);
e_approx = e00 * ones(resolution, resolution);

for i = 1:K
    for j = 1:L
        e(i,j) = (rho2*c2-rho1*c1)*(8*sqrt(Lx*Ly)/pi^2)*(1/(i*j))*cos(i*pi/2)*sin((i*pi*lx)/(2*Lx))*cos(j*pi/2)*sin((j*pi*ly)/(2*Ly));
        phi_x = eval_basis(i, x_values, Lx)';
        phi_y = eval_basis(j, y_values,Ly);
        %phi(i,j,:,:) = phi_x' * phi_y ;
        e_approx = e_approx + e(i,j)*(phi_x*phi_y);
    end
end
[X, Y] = meshgrid(x_values, y_values);
figure;
surf(X, Y, e_approx);
shading interp;  % Smooth shading
colormap jet;    % Color scheme
colorbar;        % Add color scale
xlabel('x');
ylabel('y');
zlabel('e(x,y)');
title('Approximation of e(x,y)=\rho(x,y)c(x,y)');

function e_vec = half_e(K, L_space, l)
    i_vec = 1:K;
    angle_vec = i_vec * pi / 2;
    e_vec = cos(angle_vec) .* sin(angle_vec * l/L_space) ./ angle_vec * sqrt(2* L_space);
end

function e_mat = e_mat(K,L, L_x, L_y, l_x, l_y)
    e_mat = half_e(K,L_x,l_x)' * half_e(L,L_y,l_y);
end

function phi_eval =  eval_basis(index, values, Length)
    % given values = [x1 x2 ... xn]
    % it returns [\phi_i (x1) ... \phi_i (xn)]
    % where i is the index
    % Length is the Lx or Ly used to define \phi. The point is that this
    % works for both x and y by exploiting symmetry
    phi_eval = cos(index * pi * values/Length) * sqrt(2/Length);
    if index == 0
        phi_eval = 1/sqrt(Length) * ones(length(values),1)';
    end
end