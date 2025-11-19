% Model parameters
A = 1;
% Initial condition
phi0 = 1;
% Time
t = 2;
% Force as a function
f = @(s) (s < 0).*0 + (s >= 0).*1;


% Solution
x = solution(A, phi0, f, t)



