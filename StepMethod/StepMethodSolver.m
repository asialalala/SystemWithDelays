
function [time, solution] = StepMethodSolver(tk, f_ode, Xtau, dXtau, phi)
syms y(t)

t0 = 0;
tau = max(Xtau, dXtau);

solution = [];
time = [];

phiCurrent = phi;
titer = t0;

for iter = 1:ceil(tk/tau)
    % substitution for x(t-tau)
    eq = diff(y,t) == f_ode(t, y, phiCurrent(t - Xtau), phiCurrent(t - dXtau));
    cond = y(titer) == phiCurrent(titer);
    sol = dsolve(eq, cond);
    
    phiCurrent = @(tVal) subs(sol, t, tVal); % Update phi with the latest solution

     % time span
    titer = titer + tau;
    t_start = (iter - 1) * tau;
    t_end = titer;
    time_segment = linspace(t_start, t_end, 100);
   
    % calculate
    sol_numeric = double(subs(sol, t, time_segment));
    
    % concatenate
    solution = [solution, sol_numeric];
    time = [time, time_segment]; 


end