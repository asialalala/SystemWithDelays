function [numHandle, formulaSym] = getRungeKuttaFormula()
% Derive formula for the Runge-Kutta method
%
%   Returns numHandle - numerical function, that can be used
%   for further  calculations and  formulaSym - symbolic formula.

    % 1. Symbols definition
    syms h y_n k1 k2 k3 k4 real
    
    % k1 = f(t_n, y_n)
    % k2 = f(t_n + h/2, y_n + h*k1/2)
    % k3 = f(t_n + h/2, y_n + h*k2/2)
    % k4 = f(t_n + h, y_n + h*k3)
    
    % 2. Derive symbolic formula RK4
    % formula: y_{n+1} = y_n + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    phi = (k1 + 2*k2 + 2*k3 + k4) / 6;
    formulaSym = y_n + h * phi;
    
    % 3. Define numerical function (matlabFunction)
    % variables: step h, current value y_n and four tilt
    numHandle = matlabFunction(formulaSym, 'Vars', [h, y_n, k1, k2, k3, k4]);
end