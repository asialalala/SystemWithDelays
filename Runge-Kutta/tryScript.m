clear all;
close all;

% 1. Model parameters definition
h = 0.01;
tk = 6;

A = 3;
B = 2;
C = 1;
tau = 0;
k = 2;

% 2. Problem definition
phi = exampleFunction.get_phi1();
f_ode = exampleFunction.get_f_ode1();
f_sol = exampleFunction.get_sol1();
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


