function [numHandle, formulaSym] = getImplicitAdamsFormula(k)
% Derive formula for the Implicit Adams method
% with order k.
%   k - method's order
%
%   Returns numHandle - numerical function, that can be used
%   for further  calculations and  formulaSym - symbolic formula.

    % 1. Symbols definition
    h = sym('h');
    y_n = sym('y_n');
    syms f(n);

    gamma =  getImplicitAdamsCoefficients(k); 

    % 2. Derive symbolic formula
    sum = 0;
    for j = 0:k
        sum = sum + gamma(j + 1) * backwardDifference(f, j, n+1);
    end
    formulaSym = simplify(y_n + h * sum);

    % 3. Define numerical function
    f_calls = arrayfun(@(idx) f(n+1-idx), 0:k); % [f(n), f(n-1), ..., f(n-k+1)]
    f_names = sym('f', [1 k+1]); % [f1, ..., f_k]
    cleanFormula = subs(formulaSym, f_calls, f_names);
    
    numHandle = matlabFunction(cleanFormula, 'Vars', [h, y_n, f_names]);
end