function [time, solution] = rungeKuttaSolver(k, h, tk, f_ode, tau, phi)
% Solves the given differential equation using the Runge-Kutta method.
%   h     - Time step size
%   tk    - Final simulation time
%   f_ode - Function handle defining the ODE
%   tau   - Delay from the equation
%   phi   - Initial condition function defined on the interval [-max(lags), 0]
%
%   Returns solution - solution from 0 to tk

    % 1. Runge Kutta method formula
    [rungeKuttaFunc, ~] = getRungeKuttaFormula();

    % 3. Memory grid initialization
    Ntau = round(tau/h); % Sample offset for tau (delay)
    N = round(tk/h);     % Solution number of samples  
    Ndelta = 1; % The smallest possible sample offset for derivative
    Nteta = Ntau + Ndelta;
    Nx = N + Nteta + 1;   % Number of all samples

    x = zeros(1, Nx);
    offset = tau + h;
    t = ((1:Nx) - (Nteta + 1)) * h; % time axis offset by tau
    
    % 4. Initial condition (history)
    tetaSpan = -offset: h : 0;
    x(1:Nteta+1) = phi(tetaSpan);

    % 5. Main calculation loop
    for n = Nteta+1 : Nx-1 % Start from t=0 (Ntau+1), calculate value for n+1 element
        curr_t = t(n);
        curr_x = x(n);


        % prepare directions
        k1 = f_ode(curr_t, curr_x);
        k2 = f_ode(curr_t + h/2, curr_x + h*k1/2); % uproszczenie
        k3 = f_ode(curr_t + h/2, curr_x + h*k2/2);
        k4 = f_ode(curr_t + h,   curr_x + h*k3);

        % Run: rungeKuttaFunc(h, y_n, k1,k2, ...)
        x(n+1) = rungeKuttaFunc(h, x(n), k1, k2, k3, k4);
  
    end

    % 6. Return solution from 0 to tk
    solution = x(Nteta+ 1:end);
    time = t(Nteta+1:end);
end


