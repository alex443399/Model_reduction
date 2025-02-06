rel_error_cubic_sines = load('rel errors, T0= Cubic and u= Sines.mat').error_total_rel;
rel_error_cubic_zero = load('rel errors, T0= Cubic and u= Zero.mat').error_total_rel;
rel_error_Jump_Sines = load('rel errors, T0= Jump and u= Sines.mat').error_total_rel;
rel_error_Jump_zero = load('rel errors, T0= Jump and u= Zero.mat').error_total_rel;

abs_error_cubic_sines = load('abs errors, T0= Cubic and u= Sines.mat').error_total_abs;
abs_error_cubic_zero = load('abs errors, T0= Cubic and u= Zero.mat').error_total_abs;
abs_error_Jump_Sines = load('abs errors, T0= Jump and u= Sines.mat').error_total_abs;
abs_error_Jump_zero = load('abs errors, T0= Jump and u= Zero.mat').error_total_abs;
%% Plotting relative errors
clf
semilogy([2,4,8,16,32], rel_error_cubic_zero, 'DisplayName', "T0 = Cubic, u = Zero");
hold on
semilogy([2,4,8,16,32], rel_error_cubic_sines, 'DisplayName', "T0 = Cubic, u = Sines");
hold on
semilogy([2,4,8,16,32], rel_error_Jump_zero, 'DisplayName', "T0 = Jump, u = Zero");
hold on
semilogy([2,4,8,16,32], rel_error_Jump_Sines, 'DisplayName', "T0 = Jump, u = Sines");
hold off
legend()
title("Relative errors for different experiments")
xlabel("K,L")
ylabel("e_{rel,tot}")
%% Plotting absolute errors
clf
semilogy([2,4,8,16,32], abs_error_cubic_zero, 'DisplayName', "T0 = Cubic, u = Zero");
hold on
semilogy([2,4,8,16,32], abs_error_cubic_sines, 'DisplayName', "T0 = Cubic, u = Sines");
hold on
semilogy([2,4,8,16,32], abs_error_Jump_zero, 'DisplayName', "T0 = Jump, u = Zero");
hold on
semilogy([2,4,8,16,32], abs_error_Jump_Sines, 'DisplayName', "T0 = Jump, u = Sines");
hold off
legend()
title("Absolute errors for different experiments")
xlabel("K,L")
ylabel("e_{abs,tot}")