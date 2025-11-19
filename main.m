% Model parameters
A = [1, 0; 0, 2];
% Initial condition
phi0 = [1; 1];
% Time
t = [0:2];
% Force as a function
f = @(s) (s < 0).*0 + (s >= 0).*1;
fvec = @(s) [f(s); f(s)];

% Solution
x = solution(A, phi0, fvec, t);

plot(t, x)
xlabel("time")
ylabel("x(t)")


