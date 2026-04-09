clear all;
close all;

% 1. Model parameters definition
h = 0.01;
tk = 6;

A = 3;
B = 2;
C = 1;

k = 2;

% 2. Problem definition
phi = exampleRddeFunctions.get_phi1(C);
f_ode = exampleRddeFunctions.get_f_ode1(h, A, B, C);
f_sol = exampleRddeFunctions.get_sol1(A, B, C);
tau = B;
t_span = 0:h:tk;

x_exact = f_sol(t_span); % Reference

[t_adamsS, x_adamsS] = implicitAdamsSolver(k, h, tk, f_ode, tau, phi);
error_valS = abs(x_adamsS - x_exact);

plot(t_adamsS, x_adamsS)
hold on;
grid on;
plot(t_adamsS, x_exact, "--")
xlabel("t")
ylabel("x")
legend("implicit", "exact")
title("Implicit adams method comparison with excat solution")


