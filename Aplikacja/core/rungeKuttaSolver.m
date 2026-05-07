function [time, solution] = rungeKuttaSolver(k, h, tk, f_ode, Xtau, dXtau, phi)
% Solves the given differential equation using the Runge-Kutta method.
%   h     - Time step size
%   tk    - Final simulation time
%   f_ode - Function handle defining the ODE
%   tau   - Delay from the equation
%   phi   - Initial condition function defined on the interval [-max(lags), 0]
%
%   Returns solution - solution from 0 to tk

    % 3. Memory grid initialization
    NXtau  = round(Xtau/h);
    NdXtau = round(dXtau/h);
    Ntau   = max(NXtau, NdXtau); 
    N = round(tk/h);     % Solution number of samples  

    Nteta = Ntau + 1;    % The smallest possible sample offset for derivative
    Nx = N + Nteta;  % Number of all samples

    x = zeros(1, Nx);
    t = ((1:Nx) - Nteta) * h; % time axis offset by tau
    
    % 3. Initial condition (history)
    x(1:Nteta) = phi(t(1:Nteta));
    
    % 5. Main calculation loop
    for n = Nteta : Nx-1 % Start from t=0 (Ntau+1), calculate value for n+1 element
        
        % calculate x(n+1)
        factor = getFactor(x, t, h, k, f_ode, Xtau, dXtau, n, Nteta);
        x(n+1) = x(n) + h * factor;
    end

    % 6. Return solution from 0 to tk
    solution = x(Nteta:end);
    time = t(Nteta:end);
end

function xDelayed = getXvalue(x, t_target, h, n_current, Nteta)
    
        % Map to indexes 
        idxFloat = t_target / h + Nteta;
        idxLower = floor(idxFloat);
        idxUpper = idxLower + 1;
        
        if idxUpper > n_current
            xDelayed = x(n_current);
        else
            % Linear interpolation
            alpha = idxFloat - idxLower;
            xDelayed = x(idxLower) + (x(idxUpper)-x(idxLower))*alpha;
        end
end

function dxDelayed = getDXvalue(x, t_target, h, n_current, Nteta)
    
        % Map to indexes 
        idxFloat = t_target / h + Nteta;
        idxLower = floor(idxFloat);
        idxUpper = idxLower + 1;
        
        if idxUpper > n_current
            idxUpper = n_current;
            idxLower = n_current - 1;
        end

        % Calculate the delayed derivative value
        dxDelayed = (x(idxUpper) - x(idxLower))/h;
end

function factor = getFactor(x, t, h, k, f_ode, Xtau, dXtau, n, Nteta)


    % interpolation for x()
    getD = @(t_curr) getXvalue(x, t_curr - Xtau, h, n, Nteta);
    
    % interpolation for x'()
    getXD = @(t_curr) getDXvalue(x, t_curr - dXtau, h, n, Nteta);

    switch k
        case 1
            % RK1 Euler's method
            t1 = t(n);
            k1 = f_ode(t1, x(n), getD(t1), getXD(t1));
            factor = k1;
        case 2
            % RK2 Improved Euler's method
            t1 = t(n);
            t2 = t(n) + h;
            k1 = f_ode(t1, x(n), getD(t1), getXD(t1));
            k2 = f_ode(t2, x(n) + h*k1, getD(t2), getXD(t2));
            factor = (k1 + k2) / 2;
        case 4
            % RK4
            t1 = t(n);
            t2 = t(n) + h/2;
            t3 = t(n) + h/2;
            t4 = t(n) + h;
            k1 = f_ode(t1, x(n), getD(t1), getXD(t1));
            k2 = f_ode(t2, x(n) + h*k1/2, getD(t2), getXD(t2));
            k3 = f_ode(t3, x(n) + h*k2/2, getD(t3), getXD(t3));
            k4 = f_ode(t4, x(n) + h*k3, getD(t4), getXD(t4));
            factor = (k1 + 2*k2 + 2*k3 + k4) / 6;
        case 5
            % RK 3/8
            t1 = t(n);
            t2 = t(n) + h/3;
            t3 = t(n) + 2*h/3;
            t4 = t(n) + h;
            k1 = f_ode(t1, x(n), getD(t1), getXD(t1));
            k2 = f_ode(t1 + h/3, x(n) + h*k1/3, getD(t2), getXD(t2));
            k3 = f_ode(t1 + 2*h/3, x(n) - h*k1/3 + h*k2, getD(t3), getXD(t3));
            k4 = f_ode(t1 + h, x(n) + h*k1 - h*k2 + h*k3, getD(t4), getXD(t4));
            factor = (k1 + 3*k2 + 3*k3 + k4) / 8;
        otherwise
            error('RKSolver:OrderNotSupported', ...
          'Order %d is not supported. Use 1, 2, or 4.', k);
    end
end
