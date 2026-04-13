function [time, solution] = rungeKuttaSolver(k, h, tk, f_ode, tau, phi)
% Solves the given differential equation using the Runge-Kutta method.
%   h     - Time step size
%   tk    - Final simulation time
%   f_ode - Function handle defining the ODE
%   tau   - Delay from the equation
%   phi   - Initial condition function defined on the interval [-max(lags), 0]
%
%   Returns solution - solution from 0 to tk


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
        
        % calculate x(n+1)
        factor = getFactor(x, t, h, k, f_ode, tau, n, phi);
        x(n+1) = x(n) + h * factor;
    end

    % 6. Return solution from 0 to tk
    solution = x(Nteta+ 1:end);
    time = t(Nteta+1:end);
end

function xDelayed = getValue(phi_init, x, t_target, h, n_current, Nteta)
    
    if t_target <= 0
        xDelayed = phi_init(t_target);
    else
        % Map to inexes 
        idxFloat = (t_target / h) + (Nteta + 1);
        idxLower = floor(idxFloat);
        idxUpper = idxLower + 1;
        
        if idxUpper > n_current
            xDelayed = x(n_current);
        else
            % Linear interpolation
            alpha = idxFloat - idxLower;
            xDelayed = x(idxLower) + alpha * (x(idxUpper) - idxLower);
        end
    end
end

function factor = getFactor(x, t, h, k, f_ode, tau, n, phi_init)
    Nteta = round(tau/h) + 1; 
    
    % Lokalna funkcja pomocnicza do pobierania y(t-tau)
    getD = @(t_curr) getValue(phi_init, x, t_curr - tau, h, n, Nteta);

    switch k
        case 1
            k1 = f_ode(t(n), x(n), getD(t(n)));
            factor = k1;
        case 2
            k1 = f_ode(t(n), x(n), getD(t(n)));
            k2 = f_ode(t(n) + h, x(n) + h*k1, getD(t(n) + h));
            factor = (k1 + k2) / 2;
        case 4
            % Standardowy RK4
            k1 = f_ode(t(n), x(n), getD(t(n)));
            k2 = f_ode(t(n) + h/2, x(n) + h*k1/2, getD(t(n) + h/2));
            k3 = f_ode(t(n) + h/2, x(n) + h*k2/2, getD(t(n) + h/2));
            k4 = f_ode(t(n) + h, x(n) + h*k3, getD(t(n) + h));
            factor = (k1 + 2*k2 + 2*k3 + k4) / 6;
        otherwise
            error('Obsługiwane rzędy k: 1, 2, 4.');
    end
end