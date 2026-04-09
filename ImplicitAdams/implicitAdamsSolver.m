function [time, solution] = implicitAdamsSolver(k, h, tk, f_ode, tau, phi)
% Solves the given differential equation using the Explicit Adams method.
%   k     - Method order
%   h     - Time step size
%   tk    - Final simulation time
%   f_ode - Function handle defining the ODE
%   tau   - Delay from the equation
%   phi   - Initial condition function defined on the interval [-max(lags), 0]
%
%   Returns solution - solution from 0 to tk

    % 1. Explicit Adams method formula with order k
    [exAdamsFunc, ~] = getExplicitAdamsFormula(k);

    % 2. 
    [imAdamsFunc, ~] = getImplicitAdamsFormula(k);

    % 3. Memory grid initailization
    Ntau = round(tau/h); % Sample offset fot tau (delay)
    N = round(tk/h);     % Solution number of samples  
    Ndelta = 1; % The smalest posiible sample offset for derivative
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
        
        % Check if enough history exists for order k; if not fallback to lower order (Euler), accounting for delta lag. Consider the
        % Ndelta for derivative
        k_eff = min(k, n - Ntau - Ndelta); 

        if k_eff < k 
            % STARTING  Eulera method(k = 1)
            fn = f_ode(n, x(1:n));
            x(n+1) = x(n) + h * fn;
        else
            % PREDICT ------
            % Prepare derivative vector: [f(n), f(n-1), ..., f(n-k+1)]
            f_vals = zeros(1, k); % k order -> k derivatives
            for i = 0:k-1
                idx = n - i;
                f_vals(i+1) = f_ode(idx, x(1:idx));
            end

            % Run: exAdamsFunc(h, y_n, f1, f2, ...)
            args = num2cell([h, x(n), f_vals]);
            x(n+1) = exAdamsFunc(args{:});

            % CORRECT -------
            % prepare predicted f
            f_pred = f_ode(n+1, x(1:n+1));
            % concatenated f's
            f_vals = [f_pred, f_vals];
            args = num2cell([h, x(n), f_vals]);
            x(n+1) = imAdamsFunc(args{:});
        end
  
    end

    % 6. Return solution from 0 to tk
    solution = x(Nteta+ 1:end);
    time = t(Nteta+1:end);
end


