function [time, solution] = expliciteAdamsWithoutStartSolverNA(k, h, tk, f_ode, tau, phi)
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
[adamsFunc, formulaSym] = getFormulaNA(k,h, tau);

% 2. Memory grid initailization
Ntau = round(tau/h); % Sample offset for tau (delay)
N = round(tk/h);     % Solution number of samples
Ndelta = 1; % The smalest posiible sample offset for derivative
Nteta = 2*Ntau + Ndelta + (k-1); % multiply by 2!
Nx = N + 1 + Nteta;   % Number of all samples

x = zeros(1, Nx);
offset = 2*tau + h + (k-1)*h; % multiply by 2!
t = ((1:Nx) - (Nteta + 1)) * h; % time axis offset 


% 3. Initail condition (history)
tetaSpan = -offset: h : 0;
x(1:Nteta+1) = phi(tetaSpan);

% 4. MAIN COMPUTATIONAL LOOP
for n = Nteta+1 : Nx-1 % Start od t=0 (probki Ntau+1), obliczamy wartosc dla probki n+1 probki
    
    % Check if enough history exists for order k; if not fallback to lower order (Euler), accounting for delta lag. Consider the
    % Ndelta for derivative

    % Prepare derivative vector: [f(n), f(n-1), ..., f(n-k/2+1)]
    f_vals_1 = zeros(1, ceil(k/2)); % k order -> k derivatives
    for i = 0:(round(k/2)-1)
        idx = n - i;
        f_vals_1(i+1) = f_ode(idx, x(1:idx));
    end

    % Prepare derivative vector in tau
    f_tau_vals_2 = zeros(1, floor(k/2));
    for i = 0:(floor(k/2)-1)
        idx = n - i - Ntau;
        f_tau_vals_2(i+1) = f_ode(idx, x(1:idx));
    end
    
    f_vals = [f_vals_1, f_tau_vals_2];

    % Run: adamsFunc(h, y_n, f1, f2, ...)
    args = num2cell([x(n), f_vals]);
    x(n+1) = adamsFunc(args{:});
end
% Return solution from 0 to tk
solution = x(Nteta + 1:end);
time = t(Nteta+1:end);
end