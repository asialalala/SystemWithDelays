% Model parameters
A = -1;
A1 = 1;
% Initial condition
phi = @(s) 1;
% Time
t = linspace(0, 2, 100);
% Force as a function
f = @(s) s;
% Time delay
h = 0;

% Solution
x = solution(A,  A1, h, phi, f, t);

plot(t, x)
xlabel("time")
ylabel("x(t)")


