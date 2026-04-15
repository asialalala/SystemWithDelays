function [time, solution] = explicitNddeAdamsSolver(k, h, tk, f_ode, Xtau, dXtau, phi)
% Solves the given differential equation using the Explicit Adams method.
%   k     - Method order
%   h     - Time step size
%   tk    - Final simulation time
%   f_ode - Function handle defining the ODE
%   tau   - Delay from the equation
%   phi   - Initial condition function defined on the interval [-max(lags), 0]
%
%   Returns solution - solution from 0 to tk

    % 1. Adaams method formula with order k
    [adamsFunc, formulaSym] = getExplicitAdamsFormula(k);

    % 2. Memory grid initailization
    % Sample offset fot tau (delay)
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

    % 4. MAIN COMPUTATIONAL LOOP
    for n = Nteta : Nx-1 % Start from t=0 (Ntau)
        
        % Approx derivative
        dx = approxDx(n, NdXtau, x, dXtau, h);

        % Check if enough history exists for order k; if not fallback to lower order (Euler), accounting for delta lag. Consider the
        % Ndelta for derivative
        k_eff = min(k, n - Ntau - 1); 

        if k_eff < k 
            % STARTING  Euler method(k = 1)
            fn = f_ode(t(n), x(n), x(n - NXtau), dx);
            x(n+1) = x(n) + h * fn;
        else
            % Prepare derivative vector: [f(n), f(n-1), ..., f(n-k+1)]
            f_vals = zeros(1, k); % k order -> k derivatives
            for i = 0:k-1
                idx = n - i;
                % Approx derivative
                dx = approxDx(idx, NdXtau, x, dXtau, h);
                f_vals(i+1) = f_ode(t(idx), x(idx), x(idx - NXtau), dx);
            end

            % Run: adamsFunc(h, y_n, f1, f2, ...)
            args = num2cell([h, x(n), f_vals]);
            x(n+1) = adamsFunc(args{:});
        end
    end
    % Return solution from 0 to tk
    solution = x(Nteta:end);
    time = t(Nteta:end);
end

function dx = approxDx(n, NdXtau, x, dXtau, h)
 % Approx derivative
        if dXtau == 0
            dx = 0;
        elseif n - NdXtau > 1
            dx = (x(n - NdXtau) - x(n - NdXtau - 1)) /h;
        else 
            dx = 0;
        end
end 
