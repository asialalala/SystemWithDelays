% Model parameters
A = 1;
% Initial condition
phi0 = 1;
% Time
t = [0:2];
% Force as a function
f = @(s) (s < 0).*0 + (s >= 0).*1;


% Solution
x = solution(A, phi0, f, t);

plot(t, x)
xlabel("time")
ylabel("x(t)")


