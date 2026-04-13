clear all;
close all;

% 1. Model parameters definition
h = 0.1;
tk = 10;

A = 3;
B = 2;
C = 1;
tau = B;
k = 2;

% 2. Problem definition
phi = exampleRddeFunctions.get_phi1(C);
f_ode = exampleRddeFunctions.get_f_ode1(A,B,C);
f_sol = exampleRddeFunctions.get_sol1(A, B, C);
t_span = 0:h:tk;

x_exact = f_sol(t_span); % Reference

[t_rungeS, x_rungeS] = rungeKuttaSolver(k, h, tk, f_ode, tau, phi);
error_valS = abs(x_rungeS - x_exact);

plot(t_rungeS, x_rungeS)
hold on;
grid on;
plot(t_rungeS, x_exact, "--")
xlabel("t")
ylabel("x")
legend("implicit", "exact")
title("Implicit adams method comparison with excat solution")


