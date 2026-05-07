% Lambert Function W
A = 3;
B = 2;
C_val = 1;

h = 0.1;
tk = 10;
t_span = 0:h:tk;

% 2. Problem definition
phi = exampleRddeFunctions.get_phi1(C_val);
f_sol = exampleRddeFunctions.get_sol1(A, B, C_val);

x_exact = f_sol(t_span); % Reference

% Model parameters dy(t)+alfa*y(t-T)+beta(t)=0, T>0
alfa = -A;
beta = 0;
T = B;

N = 1000; % Lambert modes used in the method

k_vec = -N:N;
s_vec = (1/T .* lambertw(k_vec, -alfa * T * exp(beta * T)) - beta);
t_vec = linspace(T, 0, 2*N + 1)' - T; % [0, N] divided on 2N elements

% Omega initialization
[S_grid, T_grid] = meshgrid(s_vec, t_vec);
Omega = exp(S_grid .* T_grid);

% Coefficients C
Phi = phi(t_vec);
C = Omega \ Phi;

[S_final, T_final] = meshgrid(s_vec, t_span);
x_matrix = exp(S_final .* T_final) .* C.'; 
x_numeric = sum(x_matrix, 2);

figure;
plot(t_span, x_exact, 'r', 'LineWidth', 2); hold on;
plot(t_span, real(x_numeric), 'b--', 'LineWidth', 1.5);
grid on;
legend('Dokładne', 'Lambert W');
title('Rozwiązanie DDE metodą Lamberta W');

