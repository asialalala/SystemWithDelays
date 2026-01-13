function [time, solution] = expliciteNddeAdamsWithoutStartSolver(k, h, tk, f_ode, tau, phi)
% Solves the given differential equation using the Explicit Adams method.
%   k     - Method order
%   h     - Time step size
%   tk    - Final simulation time
%   f_ode - Function handle defining the ODE
%   lags  - Vector of time delays (e.g., [tau, delta])
%   phi   - Initial condition function defined on the interval [-max(lags), 0]
%
%   Returns solution - solution from 0 to tk

% 1. Adaams method formula with order k
[adamsFunc, formulaSym] = getFormula(k);

% 2. Memory grid initailization
Ntau = round(tau/h); % Sample offset for tau (delay)
N = round(tk/h);     % Solution number of samples
Ndelta = 1; % The smalest posiible sample offset for derivative
Nteta = Ntau + Ndelta + (k-1);
Nx = N + 1 + Nteta;   % Number of all samples

x = zeros(1, Nx);
offset = tau + h + (k-1)*h;
t = ((1:Nx) - (Nteta + 1)) * h; % time axis offset 


% 3. Initail condition (history)
tetaSpan = -offset: h : 0;
x(1:Nteta+1) = phi(tetaSpan);

% 4. MAIN COMPUTATIONAL LOOP
for n = Nteta+1 : Nx-1 % Start od t=0 (probki Ntau+1), obliczamy wartosc dla probki n+1 probki
    
    % Check if enough history exists for order k; if not fallback to lower order (Euler), accounting for delta lag. Consider the
    % Ndelta for derivative

    % Prepare derivative vector: [f(n), f(n-1), ..., f(n-k+1)]
    f_vals = zeros(1, k); % k order -> k derivatives
    for i = 0:k-1
        idx = n - i;
        f_vals(i+1) = f_ode(idx, x(1:idx));
    end

    % Run: adamsFunc(h, y_n, f1, f2, ...)
    args = num2cell([h, x(n), f_vals]);
    x(n+1) = adamsFunc(args{:});
end
% Return solution from 0 to tk
solution = x(Nteta + 1:end);
time = t(Nteta+1:end);
end
